#pragma once

#include "PhaseFieldProblemSimple_divisionless.h"
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
    }

    virtual void initData() override
    {
      Super::initData();
      communicator_.reset(new MacroCommunicator(phaseProb_.comm(), phaseProb_.getVerletList()));
      //communicator_.reset(new MacroCommunicator(phaseProb_.comm(), phaseProb_.getVerletList(), phaseProb_.getVerletListNoOverlap()));
      interaction_.reset(new DOFVector<double>(phaseProb_.getFeSpace(0), "interaction"));
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

    // Rescaling a Phasefield from [-1,1] to [0,1] and multiplying with nematic tensor

    virtual void initTimestep(AdaptInfo* adaptInfo) override
    {
      FUNCNAME("PhaseFieldGlobal::initTimestep()");
      MSG("start (PhaseFieldGlobal::initTimestep) = %20.16e\n", mpi14::now());
      
      //phaseProb_.initTimestep(adaptInfo);

      Timer t;
      #ifdef DEBUG_PERFORMANCE
            MPI_Barrier(phaseProb_.comm());
      #endif
      MSG("time (barrier0) = %e\n", t.elapsed());

      t.reset();

      if(adaptInfo->getTimestepNumber() % 3 == 0)
      { 
        //calculating interactions only every few time steps 
        *phaseProb_.tmp_lagrange2_ << tanh(valueOf(*phaseProb_.tmp_lagrange1_) / (-std::sqrt(2.0) * eps_));
        DOFVector<double> const& phi_i = *phaseProb_.tmp_lagrange2_;
        communicator_->scatter(phi_i);
        
        // calculate interaction term
        DOFVector<double> phi_j(globalProb_.getProblem()->getFeSpace(0), "phi");
        phi_j.setCoarsenOperation(NO_OPERATION);

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
        MSG("time (gather-scatter) = %e\n", t.elapsed());
        MPI_Barrier(comm());
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
      Parameters::get("signedDist->output->write every i-th timestep",everyIthTimestep);

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
      Operator* opMuInt = new Operator(phaseProb_.getFeSpace(1));
      //scaled_In_ = (phaseProb_.current_volume_/500.0)*In_;
      addZOT(opMuInt, G(valueOf(phi))*(1.0/In_)* valueOf(*interaction_));
      // addZOT(opMuInt,(1.0/Rep_) * valueOf(*repulsion_));
      phaseProb_.getProblem()->addVectorOperator(opMuInt, 1);
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

    std::unique_ptr<DOFVector<double>> oldSolution_;

    //
    std::unique_ptr<extensions::RefinementExpression> refinement_;
    std::vector<std::pair<std::size_t,double>> contributions;

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
