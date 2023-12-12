#pragma once

#include "PhaseFieldProblemSimple_divisionless_nem_neoint.h"
#include "GlobalProblem.h"
#include "MacroCommunicator.h"
#include "ExpressionAssigner.h"
#include "InitialFunctions.h"
#include "Refinement.h"
#include "MyCouplingBaseProblem.h"


namespace AMDiS { namespace base_problems {

  class PhaseFieldGlobal
      : public CouplingBaseProblem<ProblemStat, base_problems::PhaseFieldProblem, base_problems::GlobalProblem>
  {
    using Super = CouplingBaseProblem<ProblemStat, base_problems::PhaseFieldProblem, base_problems::GlobalProblem>;

  public:
    PhaseFieldGlobal(base_problems::PhaseFieldProblem& phaseProb, base_problems::GlobalProblem& globalProb)
      : Super("main", phaseProb, globalProb)
      , phaseProb_(phaseProb)
      , globalProb_(globalProb)
    {
      Parameters::get("phase->epsilon", eps_);
      Parameters::get("phase->In", In_);
	    //Parameters::get("phase->cutTimeStep", time_threshold_);
      Parameters::get("phase->division_volume", volume_threshold_);
      Parameters::get("phase->interaction->potential", potential_);
      Parameters::get("phase->apoptosis_index", apoptosis_index_);
      Parameters::get("arh->file_prefix", arh_file_prefix_);
      Parameters::get("phase->contact_inhibiton_of_growth", contact_inhi_growth);
      Parameters::get("phase->inhibition_limit", inhibiting_strength);
      Parameters::get("phase->interaction_a", inter_a_);
      Parameters::get("phase->interaction_a_adh", a_adh);
      Parameters::get("phase->interaction_a_rep", a_rep);
    }

    virtual void initData() override
    {
      Super::initData();
      communicator_.reset(new MacroCommunicator(phaseProb_.comm(), phaseProb_.getVerletList()));
      //communicator_.reset(new MacroCommunicator(phaseProb_.comm(), phaseProb_.getVerletList(), phaseProb_.getVerletListNoOverlap()));
      interaction_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction"));
      interaction_adh_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_adh"));
      interaction_rep_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_rep"));
      interaction_Bj_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_Bj"));
      interaction_Wj_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_Wj"));
      interaction_impl_lhs_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_impl_lhs"));
      interaction_impl_rhs_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction_impl_rhs"));
	    refinement_.reset(new extensions::RefinementExpression(getProblem()->getMesh()));
      oldSolution_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "oldSolution"));

      out_ranks_.resize(comm().size());
      in_ranks_.resize(comm().size());
      for(auto out: out_ranks_){
        out.clear();
      }
      for (auto in : in_ranks_)
      {
        in.clear();
      }
      numRanks_ = int(comm().size());
      for(auto out: aliveOutStatus_){
        out.clear();
      }
      for (auto in : aliveInStatus_)
      {
        in.clear();
      }
      for(auto out: volumeOut_){
        out.clear();
      }
      for (auto in : volumeIn_)
      {
        in.clear();
      }

      if(comm().rank() < phaseProb_.numInitialCells_){
        alive = 1;
      }
    }

    // constributions to the interaction term for either attraction or repulsion
    // -----------------------------------------------------------------------------------------------

    // w'(phi_i) * B(phi_j) - repulsion
    auto wr1_B()
    {
      return [eps = eps_, potential_ = potential_](double phi_i, double phi_j) {
        phi_i = std::max(-1.0, std::min(1.0, phi_i));
        phi_j = std::max(-1.0, std::min(1.0, phi_j));
        double psi_i = (phi_i-1.0)/2;
        double psi_j = (phi_j-1.0)/2;
        double a = 1.5; //coefficient of x^4
        Parameters::get("phase->interaction_a", a);
        double b = -2.5; //coefficient of x^2
        b = -a-1;
        if(potential_==7){
          double B_j = (phi_j+1.0)/2.0;
          //double B_j = sqr(phi_j+1.0)/4.0;//modified B
          double w_der_i = 0.5*(2*b*psi_i+4*a*psi_i*sqr(psi_i));
          return phi_i > -0.99 && phi_j > -0.99
              ? B_j*w_der_i
              : 0.0;
        }
      };
    }
    // B'(phi_i) * w(phi_j) - repulsion
    auto B1_wr()
    {
      return [eps=eps_, potential_ = potential_](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double psi_i = (phi_i-1.0)/2;
          double psi_j = (phi_j-1.0)/2;
          double a = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a", a);
          double b = -2.5; //coefficient of x^2
          b = -a-1.0;
          if(potential_==7){
            double B_der_i = 0.5;
            ////double B_der_i = (phi_i+1.0)/2.0;
            double w_j = 1 + b*sqr(psi_j) + a*(sqr(sqr(psi_j)));
            return phi_i>-0.99 && phi_j > -0.99
              ?  B_der_i*w_j
              : 0.0;
          }
        };
    }

    auto phij_hat_sqr()
    {
      return [](double phi_j)
      {
        phi_j = std::max(-1.0, std::min(1.0, phi_j));
        return ((phi_j+1.0)/2.0)*((phi_j+1.0)/2.0);
      };
    }

    auto phij_doublewell()
    {
      return [](double phi_j)
      {
        phi_j = std::max(-1.0, std::min(1.0, phi_j));
        return (phi_j*phi_j-1.0)*(phi_j*phi_j-1.0);
      };
    }

    //returns (phi_i+1)/2 * (phi_j+1)/2
    auto phij_product()
    { 
      return [](double phi_i, double phi_j)
      {
        phi_j = std::max(-1.0, std::min(1.0, phi_j));
        phi_i = std::max(-1.0, std::min(1.0, phi_i));
        return (phi_i+1.0)*(phi_j+1.0)/4.0;
      };
    }

    auto Wj()
    {
      return [potential_ = potential_](double phi_j)
        {
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double psi_j = (phi_j-1.0)/2.0;
          double a = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a", a);
          double b = -2.5; //coefficient of x^2
          b = -a-1.0;
          return 1.0 + b*sqr(psi_j) + a*(sqr(sqr(psi_j)));
        };
    }

    auto Bj()
    {
      return [potential_ = potential_](double phi_j)
        {
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          return (phi_j+1.0)/2.0;
        };
    }


    // new potential repulsion implicit
    // return a_rep * sumOf((phi_j+1)**2) over all cells j
    auto neo_pot_rep_impl()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double a_r_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_rep", a_r_);
          return phi_i>-0.99 && phi_j > -0.99
              ?  a_r_*(sqr(phi_j+1))
              : 0.0;
        };
    }

    // new potential repulsion implicit
    // return a_rep * (phi_i + 1)* sumOf((phi_j+1)**2) over all cells j
    auto neo_pot_rep_expl()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double a_r_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_rep", a_r_);
          return phi_i>-0.99 && phi_j > -0.99
              ?  a_r_*(phi_i+1)*(sqr(phi_j+1))
              : 0.0;
        };
    }

    // old potential implicit rhs
    // return some long function over all cells j, check onenote page 'old interaction implicit'
    auto old_pot_impl_rhs()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double psi_i = (phi_i-1.0)/2;
          double psi_j = (phi_j-1.0)/2;
          double a = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a", a);
          double b = -2.5; //coefficient of x^2
          b = -a-1.0;
          double B_j = (phi_j+1.0)/2.0;
          double w_j = 1 + b*sqr(psi_j) + a*(sqr(sqr(psi_j)));
          double term1 = w_j/2.0;
          double term2 = B_j*(-(a+1.0)*(phi_i-1.0)/2.0 + a*(phi_i-1.0)*sqr(phi_i-1.0)/4.0);
          double term3 = -B_j*(-(a+1.0)/2.0 + (3.0*a/4.0)*sqr(phi_i-1))*phi_i;
          return phi_i>-0.99 && phi_j > -0.99
              ?  term1+term2+term3
              : 0.0;
        };
    }

            // old potential implicit lhs
    // return some long function over all cells j, check onenote page 'old interaction implicit'
    auto old_pot_impl_lhs()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double psi_i = (phi_i-1.0)/2;
          double psi_j = (phi_j-1.0)/2;
          double a = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a", a);
          double b = -2.5; //coefficient of x^2
          b = -a-1.0;
          double B_j = (phi_j+1.0)/2.0;
          double w_j = 1 + b*sqr(psi_j) + a*(sqr(sqr(psi_j)));
          double term = -B_j*(-(a+1.0)/2.0 + (3.0*a/4.0)*sqr(phi_i-1));
          return phi_i>-0.99 && phi_j > -0.99
              ?  term
              : 0.0;
        };
    }

    // new potential implicit repulsion plus adhesion lhs
    auto neo_pot_impl_lhs()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double a_r_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_rep", a_r_);
          double a_a_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_adh", a_a_);
          return phi_i>-0.99 && phi_j > -0.99
              ?  a_r_*(sqr(phi_j+1)) - 4*a_a_*(sqr(sqr(phi_j) -1.0))*(phi_i*sqr(phi_i))
              : 0.0;
        };
    }

    // new potential implicit repulsion plus adhesion rhs
    auto neo_pot_impl_rhs()
    {
      return [](double phi_i, double phi_j)
        {
          phi_i = std::max(-1.0, std::min(1.0, phi_i));
          phi_j = std::max(-1.0, std::min(1.0, phi_j));
          double a_r_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_rep", a_r_);
          double a_a_ = 1.5; //coefficient of x^4
          Parameters::get("phase->interaction_a_adh", a_a_);
          return phi_i>-0.99 && phi_j > -0.99
              ?  a_r_*(sqr(phi_j+1)) + 2*a_a_*(sqr(sqr(phi_j) -1.0))*(3.0*sqr(phi_i) - 1.0)
              : 0.0;
        };
    }

    // Rescaling a Phasefield from [-1,1] to [0,1] and multiplying with nematic tensor

    virtual void initTimestep(AdaptInfo* adaptInfo) override
    {
      FUNCNAME("PhaseFieldGlobal::initTimestep()");
      MSG("start (PhaseFieldGlobal::initTimestep) = %20.16e\n", mpi14::now());
      
      //phaseProb_.initTimestep(adaptInfo);
      
      /*
      int rank;
      MPI_Comm_rank(comm(), &rank);
      // Get the total number of nodes and name of the current node
      int comm_size;
      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Comm_size(comm(), &comm_size);
      MPI_Get_processor_name(processor_name, &name_len);

      // Print information about the current node
      printf("Rank %d running on %s\n", rank, processor_name);
      printf("Comm size%d \n", comm_size);
      */

      Timer t;
      #ifdef DEBUG_PERFORMANCE
            MPI_Barrier(phaseProb_.comm());
      #endif
      MSG("time (barrier0) = %e\n", t.elapsed());

      t.reset();
      //calculating interactions only every few time steps 
      if(adaptInfo->getTimestepNumber() % 1 == 0)
      { 
        *phaseProb_.tmp_lagrange2_ << tanh(valueOf(*phaseProb_.tmp_lagrange1_) / (-std::sqrt(2.0) * eps_));
        DOFVector<double> const& phi_i = *phaseProb_.tmp_lagrange2_;
        communicator_->scatter(phi_i);
        // calculate interaction term
        DOFVector<double> phi_j(globalProb_.getProblem()->getFeSpace(0), "phi");
        phi_j.setCoarsenOperation(NO_OPERATION);

        if(potential_ == 7)
        {
          // Gather interaction terms
          interaction_->set(0.0);
          auto& interaction = *interaction_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction, b1_w=this->B1_wr(), w1_b=this->wr1_B()](int macroIndex, DOFVector<double> const& phi_j)
            {
              return transfer_mm(phi_i, phi_j, interaction, macroIndex, [b1_w,w1_b](double p_i, double p_j)
              {
                return b1_w(p_i,p_j) + w1_b(p_i,p_j);
              });
            });         

        }
        if((potential_ == 8) || (potential_ == 10))
        {
          // Gather interaction terms
          //store product of rescaled phi_i and phi_j
          interaction_->set(0.0);
          auto& interaction = *interaction_;
          //store double well sums of all cells j
          interaction_adh_->set(0.0);
          auto& interaction_adh = *interaction_adh_;
          //store phiJ_hat_sqr sum of all cells j
          interaction_rep_->set(0.0);
          auto& interaction_rep = *interaction_rep_;


          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_adh, &interaction_rep, &interaction, pj_sqr=this->phij_hat_sqr(), pj_dw=this->phij_doublewell(), phij_prod=this->phij_product()](int macroIndex, DOFVector<double> const& phi_j)
            {
              transfer_mm_dual_unary_part2(phi_j, interaction_adh, interaction_rep, macroIndex, [pj_dw](double p_j)
              {
                return pj_dw(p_j);
              }, [pj_sqr](double p_j)
              {
                return pj_sqr(p_j);
              });
              return transfer_mm(phi_i, phi_j, interaction, macroIndex, [phij_prod](double p_i, double p_j)
              {
                return phij_prod(p_i, p_j);
              });
            });  
              
        }

        if (potential_ == 9){

          // Gather interaction terms
          //store product of rescaled phi_i and phi_j
          interaction_->set(0.0);
          auto& interaction = *interaction_;
          //store (phij+1)/2 sums of all cells j
          interaction_Bj_->set(0.0);
          auto& interaction_Bj = *interaction_Bj_;
          //store sum of w(j) for all cells j
          interaction_Wj_->set(0.0);
          auto& interaction_Wj = *interaction_Wj_;


          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_Wj, &interaction_Bj, &interaction, wj_=this->Wj(), bj_=this->Bj(), phij_prod=this->phij_product()](int macroIndex, DOFVector<double> const& phi_j)
          {
            transfer_mm_dual_unary (phi_i, phi_j, interaction_Wj, interaction_Bj, macroIndex, [wj_](double p_j)
            {
              return wj_(p_j);
            }, [bj_](double p_j)
            {
              return bj_(p_j);
            });

            return transfer_mm(phi_i, phi_j, interaction, macroIndex, [phij_prod](double p_i, double p_j)
            {
              return phij_prod(p_i, p_j);
            });
          });  
          
        }

        if(potential_ == 11) //new potential repulsion only implicit, all calculations inside the transfermm and also selection limit for phi_i and phi_j at -0.99
        {
          // Gather interaction terms
          interaction_->set(0.0);
          auto& interaction = *interaction_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction, neo_pot_rep_impl_=this->neo_pot_rep_impl()](int macroIndex, DOFVector<double> const& phi_j)
            {
              return transfer_mm(phi_i, phi_j, interaction, macroIndex, [neo_pot_rep_impl_](double p_i, double p_j)
              {
                return neo_pot_rep_impl_(p_i,p_j);
              });
            });         

        }

        if(potential_ == 12) //new potential repulsion only explicit, all calculations inside the transfermm and also selection limit for phi_i and phi_j at -0.99
        {
          // Gather interaction terms
          interaction_->set(0.0);
          auto& interaction = *interaction_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction, neo_pot_rep_expl_=this->neo_pot_rep_expl()](int macroIndex, DOFVector<double> const& phi_j)
            {
              return transfer_mm(phi_i, phi_j, interaction, macroIndex, [neo_pot_rep_expl_](double p_i, double p_j)
              {
                return neo_pot_rep_expl_(p_i,p_j);
              });
            });         

        }

        if(potential_ == 13) //old potential implicit everything in transfer_mm
        {
          // Gather interaction terms
          interaction_impl_lhs_->set(0.0);
          interaction_impl_rhs_->set(0.0);
          auto& interaction_impl_lhs = *interaction_impl_lhs_;
          auto& interaction_impl_rhs = *interaction_impl_rhs_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_impl_lhs, &interaction_impl_rhs, old_pot_impl_lhs_=this->old_pot_impl_lhs(), old_pot_impl_rhs_=this->old_pot_impl_rhs()](int macroIndex, DOFVector<double> const& phi_j)
            {
              auto temp_something_ = transfer_mm(phi_i, phi_j, interaction_impl_rhs, macroIndex, [old_pot_impl_rhs_](double p_i, double p_j)
              {
                return old_pot_impl_rhs_(p_i,p_j);
              });

              return transfer_mm(phi_i, phi_j, interaction_impl_lhs, macroIndex, [old_pot_impl_lhs_](double p_i, double p_j)
              {
                return old_pot_impl_lhs_(p_i,p_j);
              });
            });             

        }

        if(potential_ == 14) //new potential implicit everything in transfer_mm
        {
          // Gather interaction terms
          interaction_impl_lhs_->set(0.0);
          interaction_impl_rhs_->set(0.0);
          auto& interaction_impl_lhs = *interaction_impl_lhs_;
          auto& interaction_impl_rhs = *interaction_impl_rhs_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_impl_lhs, &interaction_impl_rhs, neo_pot_impl_lhs_=this->neo_pot_impl_lhs(), neo_pot_impl_rhs_=this->neo_pot_impl_rhs()](int macroIndex, DOFVector<double> const& phi_j)
            {
              auto temp_something_ = transfer_mm(phi_i, phi_j, interaction_impl_rhs, macroIndex, [neo_pot_impl_rhs_](double p_i, double p_j)
              {
                return neo_pot_impl_rhs_(p_i,p_j);
              });

              return transfer_mm(phi_i, phi_j, interaction_impl_lhs, macroIndex, [neo_pot_impl_lhs_](double p_i, double p_j)
              {
                return neo_pot_impl_lhs_(p_i,p_j);
              });
            });             

        }

        if(potential_ == 15) //old potential implicit everything in transfer_mm
        {
          // Gather interaction terms
          interaction_impl_lhs_->set(0.0);
          interaction_impl_rhs_->set(0.0);
          auto& interaction_impl_lhs = *interaction_impl_lhs_;
          auto& interaction_impl_rhs = *interaction_impl_rhs_;
          
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_impl_lhs, &interaction_impl_rhs, old_pot_impl_lhs_=this->old_pot_impl_lhs(), old_pot_impl_rhs_=this->old_pot_impl_rhs()](int macroIndex, DOFVector<double> const& phi_j)
            {
              return transfer_mm_dual_binary_part2(phi_i, phi_j, interaction_impl_rhs, interaction_impl_lhs, macroIndex, [old_pot_impl_rhs_](double p_i, double p_j)
              {
                return old_pot_impl_rhs_(p_i,p_j);
              }, [old_pot_impl_lhs_](double p_i, double p_j)
              {
                return old_pot_impl_lhs_(p_i,p_j);
              });
            });             

        }
        if(potential_ == 16) //new potential implicit everything in transfer_mm
        {
          // Gather interaction terms
          interaction_impl_lhs_->set(0.0);
          interaction_impl_rhs_->set(0.0);
          auto& interaction_impl_lhs = *interaction_impl_lhs_;
          auto& interaction_impl_rhs = *interaction_impl_rhs_;
          // Store contributions by repulsion term for each other cell
          // to clarify which one is actually a neighbour
          //std::vector<std::pair<std::size_t,double>> contributions;
          contributions.clear();

          communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction_impl_lhs, &interaction_impl_rhs, neo_pot_impl_lhs_=this->neo_pot_impl_lhs(), neo_pot_impl_rhs_=this->neo_pot_impl_rhs()](int macroIndex, DOFVector<double> const& phi_j)
            {
              return transfer_mm_dual_binary_part2(phi_i, phi_j, interaction_impl_rhs, interaction_impl_lhs, macroIndex, 
              [neo_pot_impl_rhs_](double p_i, double p_j)
              {
                return neo_pot_impl_rhs_(p_i,p_j);
              }, 
              [neo_pot_impl_lhs_](double p_i, double p_j)
              {
                return neo_pot_impl_lhs_(p_i,p_j);
              });
            });   
        }
        MSG("time (gather-scatter) = %e\n", t.elapsed());
        //MPI_Barrier(comm());
      }

      // ----- Neighbour and Interaction output from here ------------
      std::string directory_ = ".";
      std::string postfix_ = "";
      Parameters::get("output", directory_);
      Parameters::get("postfix", postfix_);
      postfix_ += "_p" + std::to_string(comm().rank());

      int writeNeighbours = 0;
      Parameters::get("main->write neighbors",writeNeighbours);

      
      int everyIthTimestep = 1;
      //Parameters::get("phase->space->output->write every i-th timestep:",everyIthTimestep);

      if (writeNeighbours && (adaptInfo->getTimestepNumber() % everyIthTimestep == 0)){
        std::string filename = directory_ + "/neighbours" + postfix_ + ".dat";
        std::ofstream neighbourFile;
        neighbourFile.open(filename,std::ios::app);
        for (std::size_t i = 0; i < contributions.size(); i++) 
          if (contributions[i].second > 0.001) 
            neighbourFile << contributions[i].first << " ";
        neighbourFile <<  "\n";
        neighbourFile.close();  
      }

      int writeInteraction = 0;
      Parameters::get("main->write interaction term",writeInteraction);
      if (writeInteraction){
        std::string filename = directory_ + "/interaction" + postfix_ + ".dat";
        std::ofstream interactionFile;
        interactionFile.open(filename,std::ios::app);
        interactionFile << "-------------New timestep number: " << adaptInfo->getTime() << "-------------------" << "\n";
        for (std::size_t i = 0; i < contributions.size(); i++){
          interactionFile << "Cell: " << contributions[i].first << "\n";
          double mu_integral = contributions[i].second/In_;
          interactionFile << "Interaction in l1 norm: " << mu_integral << "\n";
        }
        interactionFile.close();
      } 

      int writeARH = 1;
      int writeARHEveryTimeStep = 10000;
      Parameters::get("sendDOF->space->output->write every i-th timestep",writeARHEveryTimeStep);
      Parameters::get("sendDOF->space->output->write arh",writeARH);
      if(writeARH && (adaptInfo->getTimestepNumber() % writeARHEveryTimeStep == 0) && (adaptInfo->getTimestepNumber()>0))
      {
        int time = adaptInfo->getTimestepNumber();
        char time_step[10];
        sprintf(time_step, "%06d", time);
        DOFVector<double> *phasefieldDOF; // = phaseProb_.getProblem()->getSolution(0);
        phasefieldDOF = new DOFVector<double>(phaseProb_.getFeSpace(0), "phasefield");
        *phasefieldDOF << valueOf(phaseProb_.getProblem()->getSolution(0));//
        std::ofstream arhFile;
        std::cerr << "Printing arh";
        std::string arh_filename{directory_ + "/t_" + std::string(time_step) + postfix_ + ".arh"};
        //std::string arh_filename = directory_ + "/arh/t" + std::string(time_step) + postfix_ + ".arh";
        arhFile.open(arh_filename,std::ios::app);
        AMDiS::io::Arh3Writer::writeFile(phasefieldDOF, arh_filename);
        arhFile.close();
      }

      if(contact_inhi_growth == 1){
        double total_contributions;
        phaseProb_.neighbours = 0;
        for (std::size_t i = 0; i < contributions.size(); i++){
          total_contributions += contributions[i].second;
          //if(abs(contributions[i].second) > 0.001){//non zero contributions means that a cell has a neighbour
          //  phaseProb_.neighbours += 1;
          //}
        }
        total_contributions = total_contributions/In_;        
        phaseProb_.total_interactions = total_contributions;
        phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_*(std::max(0.0, std::min(1.0 - ((total_contributions > 0) - (total_contributions < 0))*std::pow((total_contributions/inhibiting_strength),2),1.0)));
        //phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_*(std::max(0.0, std::min(1.0 - (total_contributions/std::abs(total_contributions))*std::pow((total_contributions/inhibiting_strength),2),1.0)));
        //if(phaseProb_.current_volume_ < volume_threshold_*0.5){   //if the cell is small, there is no contact inhibition
        //  phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_; 
        //}
      }     

      phaseProb_.initTimestep(adaptInfo);
      globalProb_.initTimestep(adaptInfo); //why did I remove this line before?
    }

      //Stabilization
    template <class Phi>
    auto G(Phi &&p)
    { 
      return sqrt(4.0/9.0*pow<2>(std::forward<Phi>(p)-1.0)*pow<2>(std::forward<Phi>(p)+1.0));
      //return sqrt(4.0/9.0*pow<2>(std::forward<Phi>(p)-1.0)*pow<2>(std::forward<Phi>(p)+1.0)+eps_*eps_*1.0e-1);
    }

    virtual void fillCouplingOperators() override
    {
      auto *phi = phaseProb_.getProblem()->getSolution(0);
      auto phi_restrict = max(-1.0, min(1.0, valueOf(phi)));

      if (potential_ == 7) { // master thesis potential explicit
        Operator* opMuInt = new Operator(phaseProb_.getFeSpace(1));
        //scaled_In_ = (phaseProb_.current_volume_/500.0)*In_;
        addZOT(opMuInt, G(valueOf(phi))*(1.0/In_)* valueOf(*interaction_));
        // addZOT(opMuInt,(1.0/Rep_) * valueOf(*repulsion_));
        phaseProb_.getProblem()->addVectorOperator(opMuInt, 1);
      }
      if(potential_ == 8) {
        //adding operators for implicit repulsion
        Operator* opIntRep1 = new Operator(phaseProb_.getFeSpace(1), phaseProb_.getFeSpace(0));
        addZOT(opIntRep1, G(phi_restrict)*(1.0/In_)*(-a_rep)*valueOf(*interaction_rep_));
        phaseProb_.getProblem()->addMatrixOperator(opIntRep1, 1, 0);
        Operator* opIntRep2 = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntRep2, G(phi_restrict)*(1.0/In_)*(a_rep)*valueOf(*interaction_rep_));
        phaseProb_.getProblem()->addVectorOperator(opIntRep2, 1);

        //adding operators for implicit adhesion
        Operator* opIntAdh1 = new Operator(phaseProb_.getFeSpace(1), phaseProb_.getFeSpace(0));
        addZOT(opIntAdh1, G(phi_restrict)*(1.0/In_)*(2.0*a_adh)*valueOf(*interaction_adh_)*(3.0*phi_restrict*phi_restrict-1.0));
        phaseProb_.getProblem()->addMatrixOperator(opIntAdh1, 1, 0);
        Operator* opIntAdh2 = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntAdh2, G(phi_restrict)*(1.0/In_)*(4.0*a_adh)*valueOf(*interaction_adh_)*(phi_restrict*phi_restrict*phi_restrict));
        phaseProb_.getProblem()->addVectorOperator(opIntAdh2, 1);
      }
      if(potential_ == 9){//master thesis potential implicit
        Operator* opIntmuderLHS = new Operator(phaseProb_.getFeSpace(1), phaseProb_.getFeSpace(0));
        addZOT(opIntmuderLHS, G(phi_restrict)*(-1.0/In_)*(-0.5*(inter_a_+1) + (0.75*inter_a_)*(phi_restrict-1.0)*(phi_restrict-1.0))*valueOf(*interaction_Bj_));
        phaseProb_.getProblem()->addMatrixOperator(opIntmuderLHS, 1, 0);
        Operator* opIntmuderRHS = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntmuderRHS, G(phi_restrict)*(-1.0/In_)*(-0.5*(inter_a_+1) + (0.75*inter_a_)*(phi_restrict-1.0)*(phi_restrict-1.0))*valueOf(*interaction_Bj_)*phi_restrict);
        phaseProb_.getProblem()->addVectorOperator(opIntmuderRHS, 1); //also on right side

        Operator* opIntmu = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntmu, G(phi_restrict)*(1.0/In_)*((0.5*valueOf(*interaction_Wj_)) + ((valueOf(*interaction_Bj_))*((1-phi_restrict)*(0.5*(inter_a_+1))+(inter_a_/4.0)*(phi_restrict-1.0)*(phi_restrict-1.0)*(phi_restrict-1.0)))));
        phaseProb_.getProblem()->addVectorOperator(opIntmu, 1);
      }

      if (potential_ == 10){ //explicit repulsion and implicit adhesion similar to potential 7
        // operators for explictit repulsion
        Operator* opIntRep2 = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntRep2, G(phi_restrict)*(1.0/In_)*(a_rep)*valueOf(*interaction_rep_)*(phi_restrict+1));
        phaseProb_.getProblem()->addVectorOperator(opIntRep2, 1);

        //adding operators for implicit adhesion
        //Operator* opIntAdh1 = new Operator(phaseProb_.getFeSpace(1));
        //addZOT(opIntAdh1, G(phi_restrict)*(1.0/In_)*(2*a_adh)*valueOf(*interaction_adh_)*(3*phi_restrict*phi_restrict-1));
        //phaseProb_.getProblem()->addMatrixOperator(opIntAdh1, 1, 0);
        //Operator* opIntAdh2 = new Operator(phaseProb_.getFeSpace(1));
        //addZOT(opIntAdh2, G(phi_restrict)*(1.0/In_)*(4*a_adh)*valueOf(*interaction_adh_)*(phi_restrict*phi_restrict*phi_restrict));
        //phaseProb_.getProblem()->addVectorOperator(opIntAdh2, 1);
      }
      if (potential_ == 11) {//new potential repulsion only implicit, all calculations inside the transfermm and also selection limit for phi_i and phi_j at -0.99
        //adding operators for implicit repulsion
        Operator* opIntRep1 = new Operator(phaseProb_.getFeSpace(1), phaseProb_.getFeSpace(0));
        addZOT(opIntRep1, G(phi_restrict)*(-1.0/In_)*valueOf(*interaction_));
        phaseProb_.getProblem()->addMatrixOperator(opIntRep1, 1, 0);
        Operator* opIntRep2 = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntRep2, G(phi_restrict)*(1.0/In_)*valueOf(*interaction_));
        phaseProb_.getProblem()->addVectorOperator(opIntRep2, 1);
      }
      if (potential_ == 12) {//new potential repulsion only explicit, all calculations inside the transfermm and also selection limit for phi_i and phi_j at -0.99
        //adding operators for explicit repulsion
        Operator* opIntRep2 = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opIntRep2, G(phi_restrict)*(1.0/In_)*valueOf(*interaction_));
        phaseProb_.getProblem()->addVectorOperator(opIntRep2, 1);
      }
      if ((potential_ == 13) || (potential_ == 14) || (potential_ == 15) || (potential_ == 16)) { // master thesis potential implicit
        Operator* opMuIntrhs = new Operator(phaseProb_.getFeSpace(1));
        addZOT(opMuIntrhs, G(valueOf(phi))*(1.0/In_)*valueOf(*interaction_impl_rhs_));
        phaseProb_.getProblem()->addVectorOperator(opMuIntrhs, 1);
        Operator* opMuIntlhs = new Operator(phaseProb_.getFeSpace(1), phaseProb_.getFeSpace(0));
        addZOT(opMuIntlhs, G(valueOf(phi))*(1.0/In_)*valueOf(*interaction_impl_lhs_));
        phaseProb_.getProblem()->addMatrixOperator(opMuIntlhs, 1, 0);
      }

    }
    
    DOFVector<double>* getPhi()
    {
      return phaseProb_.getProblem()->getSolution(0);
    }
    
    WorldVector<double> getPosition() const
    {
      return phaseProb_.getPosition();  
    }    
    
    // return center of cell
    WorldVector<double> getCenter() const
    {
      return phaseProb_.getCenter();
    }
    
    // return upper bound of radius of cell
    double getRadius() const
    {
      return phaseProb_.getRadius();
    }
    
    mpi14::Communicator const& comm() const { return phaseProb_.comm(); }
    mpi14::Communicator& comm() { return phaseProb_.comm(); }
    
  private:
    base_problems::PhaseFieldProblem& phaseProb_;
    base_problems::GlobalProblem& globalProb_;

    std::unique_ptr<MacroCommunicator> communicator_;
    std::vector<DOFVector<double>*> phis_;

    std::unique_ptr<DOFVector<double>> interaction_; 
    std::unique_ptr<DOFVector<double>> interaction_adh_; 
    std::unique_ptr<DOFVector<double>> interaction_rep_; 
    std::unique_ptr<DOFVector<double>> interaction_Bj_; 
    std::unique_ptr<DOFVector<double>> interaction_Wj_; 
    std::unique_ptr<DOFVector<double>> interaction_impl_rhs_; 
    std::unique_ptr<DOFVector<double>> interaction_impl_lhs_; 

    double a_rep = 10.0;
    double a_adh = 0.0;

    std::unique_ptr<DOFVector<double>> oldSolution_;

    //
    std::unique_ptr<extensions::RefinementExpression> refinement_;
    std::vector<std::pair<std::size_t,double>> contributions;
    std::vector<std::pair<std::size_t,double>> Bi_Bj_contributions;

    double eps_ = 0.1;

    double In_ = 1.0;   // Interaction scaling 
    double scaled_In_ = 1.0; //scaling interaction with volume
    int potential_ = 0; // Interaction potential
	  int time_threshold_ = 10;
    int time_counter_ = 0;
    double volume_threshold_ = 491.0;
    double thetaCut;

    int contact_inhi_growth = 1;
    double inhibiting_strength = 3000.0;

    bool just_divided_ = false;
    bool just_born_ = false;
    std::vector<std::vector<int>> out_ranks_; //block chain like structure for cells which should receive DOF from a given rank
    std::vector<std::vector<int>> in_ranks_;
    std::vector<bool> out_status;
    std::vector<bool> in_status;
    int numRanks_;
    int status = 0;
    int apoptosis_index_ = 10; //generation at which cell shrinks and dies
    std::string arh_file_prefix_;

    int mother_index = 0;
    int daughter_index = 1;
    int alive = 0;

    std::vector<std::vector<int>> aliveOutStatus_;
    std::vector<std::vector<int>> aliveInStatus_;
    std::vector<std::vector<float>> volumeOut_;
    std::vector<std::vector<float>> volumeIn_;
    
    double inter_a_ = 0.001; //coefficent 'a' of interaction potential
  };
} } // end namespaces
