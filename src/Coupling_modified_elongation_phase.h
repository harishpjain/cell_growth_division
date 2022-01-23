#pragma once

#include "PhaseFieldProblemSimple_modified_elongation.h"
#include "GlobalProblem.h"
#include "MacroCommunicator.h"
#include "ExpressionAssigner.h"
#include "InitialFunctions.h"
#include "Refinement.h"
#include "MyCouplingBaseProblem.h"
//#include "FileWriter2.h"

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
        double a = 1.0;//1.5; //coefficient of x^4
        double b = -2.0;//-2.5; //coefficient of x^2
        if(potential_==7){
          double B_j = (phi_j+1.0)/2.0;
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
          double a = 1.0;//1.; //coefficient of x^4
          double b = -2.0;//-2.5; //coefficient of x^2
          if(potential_==7){
            double B_der_i = 0.5;
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
      
      int rank_ = comm().rank();
      ///time_counter_++;
      
      out_ranks_.clear();
      out_ranks_.resize(comm().size());
      in_ranks_.clear();
      in_ranks_.resize(comm().size());
      
      aliveInStatus_.clear();
      aliveInStatus_.resize(comm().size());
      aliveOutStatus_.clear();
      aliveOutStatus_.resize(comm().size());

      volumeIn_.clear();
      volumeIn_.resize(comm().size());
      volumeOut_.clear();
      volumeOut_.resize(comm().size());
      phaseProb_.volumeList_.clear();
      phaseProb_.volumeList_.resize(comm().size(), 0.0);

      phaseProb_.aliveCells_.clear();
      phaseProb_.waitingCells_.clear();
      phaseProb_.aliveCells_.resize(comm().size());
      phaseProb_.waitingCells_.resize(comm().size(),comm().size()-1);

      aliveOutStatus_[rank_].push_back(alive);
      volumeOut_[rank_].push_back(phaseProb_.volume_);

      int dummy  = 1;
      

      //making sure all cores know which other cores have cells. 
      communicator_->sendRanks(aliveOutStatus_, numRanks_, rank_); //the send and recv rank functions can also send/recv aliveStatus
      communicator_->recvRanks(aliveInStatus_, numRanks_, rank_);

      //making sure all cores have volume of cells of all other cores.
      communicator_->sendVolume(volumeOut_, numRanks_, rank_); 
      communicator_->recvVolume(volumeIn_, numRanks_, rank_);
      
      aliveInStatus_[rank_][0] = alive; //this is necessary for some reason becuase self status is not well communicated
      volumeIn_[rank_][0] = phaseProb_.volume_;

      int idx1 = 0;
      int idx2 = 0;
      
      for(int r=0; r<aliveInStatus_.size(); ++r){
        
        if(aliveInStatus_[r][0] == 1 || r < phaseProb_.numInitialCells_){
          phaseProb_.aliveCells_[idx1] = r;
          phaseProb_.volumeList_[idx1] = volumeIn_[r][0];
          idx1++;
        }else if(aliveInStatus_[r][0] == 0){
          phaseProb_.waitingCells_[idx2] = r;
          idx2++;
        }
      }

      //sorting waitingCells based on decreasing volume: NEED MORE EFFICIENT WAY TO DO THIS
      
      std::vector<std::size_t> index_vec;

      for (std::size_t i = 0; i != phaseProb_.waitingCells_.size(); ++i) { index_vec.push_back(i); }
      std::sort(
          index_vec.begin(), index_vec.end(),
          [&](std::size_t a, std::size_t b) { 
              return phaseProb_.volumeList_[a] > phaseProb_.volumeList_[b];
             });

      for( int i = 0; i < phaseProb_.waitingCells_.size() - 1; ++i )
      { 
          // while the element i is not yet in place 
          while( i != index_vec[i] )
          {
              // swap it with the element at its final place
              int alt = index_vec[i];
              std::swap( phaseProb_.waitingCells_[i], phaseProb_.waitingCells_[alt] );
              std::swap( index_vec[i], index_vec[alt] );
          }
      }
      

      for(int i = 0; i < phaseProb_.waitingCells_.size(); i++){
        if(rank_ == phaseProb_.waitingCells_[i]){
          phaseProb_.parent_ = phaseProb_.aliveCells_[i];
        }
      }

      //std::cerr << phaseProb_.aliveCells_ << "\n";
      MPI_Barrier(comm());
      if(phaseProb_.volume_ >= volume_threshold_){
        /*
        First make a rectangle cut
        Then save one half in a file, and change current phase field to half
        */

        just_divided_ = true;  //indication that the cell has just divided. Useful for signedDist Calc later
        phaseProb_.div_refine_ = 10;
        std::cerr << "Dividing Cell "<< rank_ << "\n";
        //refinement_->refine(30, function_(indicator3(getMesh()->getName() + "->phase"/*, 5 * eps_*/),
        //                                      valueOf(*phaseProb_.getProblem()->getSolution(0))));
        refinement_->refine(5, function_(indicator4(getMesh()->getName() + "->phase", 5 * eps_),
                                              valueOf(*phaseProb_.tmp_lagrange1_)));
        for (int _rep = 0; _rep < 5; ++_rep)
        {
          phaseProb_.calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(3, function_(indicator4(getMesh()->getName() + "->phase", 7 * eps_),
                                          valueOf(*phaseProb_.tmp_lagrange1_)));
        }
        //thetaCut = rand() % 180 + 1;
        thetaCut = phaseProb_.major_axis_angle_;
        WorldVector<double> r_center_ = integrate(X() * ((valueOf(*phaseProb_.tmp_lagrange2_) + 1.0) / 2.0)) / integrate((valueOf(*phaseProb_.tmp_lagrange2_) + 1.0) / 2.0); //Center of Mass
        WorldVector<double> domainDimension_;
        Parameters::get(getMesh()->getName() + "->dimension", domainDimension_);
        *phaseProb_.getProblem()->getSolution(0) << ((0.5 * (valueOf(phaseProb_.getProblem()->getSolution(0)) + 1.0)) * function_(rectangleCut(eps_, thetaCut, r_center_, domainDimension_[0] * 2), X()) * 2.0) - 1.0;
        //*phaseProb_.getProblem()->getSolution(0) << max(-1.0, min(1.0, valueOf(*phaseProb_.getProblem()->getSolution(0))));
        *phaseProb_.getProblem()->getSolution(0) << tanh(valueOf(*phaseProb_.getProblem()->getSolution(0)));
        for (int _rep = 0; _rep < 5; ++_rep)
        {
          phaseProb_.calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          //refinement_->refine(10, function_(indicator3(getMesh()->getName() + "->phase", 5 * eps_),
          //                                valueOf(*phaseProb_.getProblem()->getSolution(0))));
          refinement_->refine(10, function_(indicator4(getMesh()->getName() + "->phase", 5 * eps_),
                                              valueOf(*phaseProb_.tmp_lagrange1_)));
                                              
          *phaseProb_.getProblem()->getSolution(0) << (((0.5 * (valueOf(phaseProb_.getProblem()->getSolution(0)) + 1.0)) * function_(rectangleCut(eps_, thetaCut, r_center_, domainDimension_[0
                  ] * 2), X()) * 2.0) - 1.0);
          //*phaseProb_.getProblem()->getSolution(0) << max(-1.0, min(1.0, valueOf(*phaseProb_.getProblem()->getSolution(0))));
          *phaseProb_.getProblem()->getSolution(0) << tanh(valueOf(*phaseProb_.getProblem()->getSolution(0)));
        }

        
        DOFVector<double> *otherHalf; // = phaseProb_.getProblem()->getSolution(0);
        otherHalf = new DOFVector<double>(phaseProb_.getFeSpace(0), "other pahse");
        //*otherHalf << valueOf(phaseProb_.getProblem()->getSolution(0));//
        *otherHalf << ((0.5 * (0.5 * valueOf(phaseProb_.getProblem()->getSolution(0)) + 0.5)) * (function_(halfCutTwo(eps_, thetaCut, r_center_), X())) - 0.5) / 0.5;
        DOFVector<double> const &otherHalfRef = *otherHalf;

        //experimental
        //non iterator method
        std::cerr << phaseProb_.aliveCells_ << "\n";
        std::cerr << phaseProb_.waitingCells_ << "\n";
        std::cerr << phaseProb_.volumeList_ << "\n";
        
        for(int i = 0; i < phaseProb_.aliveCells_.size(); i++){
          if(rank_ == phaseProb_.aliveCells_[i]){
            mother_index = i;
            break;
          }
        }
        
        //auto indexItr = std::find(phaseProb_.aliveCells_.begin(), phaseProb_.aliveCells_.end(), rank_); //finds index of mother cell in alive cells
        //auto mother_index = std::distance(phaseProb_.aliveCells_.begin(), indexItr);
        int daughter_rank = phaseProb_.waitingCells_[mother_index];
        std::cerr << "\n I am " << rank_;
        std::cerr << "and my daughter is " << daughter_rank;
        std::cerr << "      *****\n";
        char from_r[10];
        sprintf(from_r, "%03d", rank_);
        char to_r[10];
        //sprintf(to_r, "%03d", phaseProb_.rankTree_[phaseProb_.nextChildIndex_]);
        //int child = phaseProb_.rankTree_[phaseProb_.nextChildIndex_];
        //out_ranks_[rank_].push_back(child);
        sprintf(to_r, "%03d", daughter_rank);//phaseProb_.waitingCells_[mother_index]);
        out_ranks_[rank_].push_back(daughter_rank);//phaseProb_.waitingCells_[mother_index]);

        ///std::string poststr = str(from_r) + "_to_rank_" + str(to_r);
        phaseProb_.nextChildIndex_++;
        if(phaseProb_.nextChildIndex_== apoptosis_index_){
          phaseProb_.growth_factor_ = -phaseProb_.growth_factor_;
        }
        std::string arh_filename{arh_file_prefix_ + std::string(from_r) + "_to_rank_" + std::string(to_r) + ".arh"};
        AMDiS::io::Arh3Writer::writeFile(otherHalf, arh_filename);//"output/data/rank_" + std::string(from_r) + "_to_rank_" + std::string(to_r) + ".arh");
        
        *phaseProb_.getProblem()->getSolution(0) << ((0.5 * (0.5 * valueOf(phaseProb_.getProblem()->getSolution(0)) + 0.5)) * (function_(halfCutOne(eps_, thetaCut, r_center_), X())) - 0.5) / 0.5; //Cutting the cell
        delete otherHalf;
      }

      //MPI_Barrier(comm());
      communicator_->sendRanks(out_ranks_, numRanks_, rank_);
      communicator_->recvRanks(in_ranks_, numRanks_, rank_);
      //MPI_Barrier(comm());
    
      status = 0;
      
      for (auto i : in_ranks_)
      {
        status += static_cast<int>(i.size());
      }
      phaseProb_.phaseStatus = status;
      if (status > 0 && just_born_ == false)
      {
        for (int in_rank : in_ranks_[phaseProb_.parent_])
        {
          if (rank_ == in_rank)
          {
            just_born_ = true;
            phaseProb_.div_refine_ = 10;
            phaseProb_.alive_ = 1;
            alive = 1;
            DOFVector<double> phi_temp(phaseProb_.getProblem()->getFeSpace(0), "phi"); //empty phase to temporarily store DOF
            //phi_temp.setCoarsenOperation(NO_OPERATION);
            //DOFVector<double> phi_temp(globalProb_.getProblem()->getFeSpace(0), "phi");
            //phi_temp = *phaseProb_.getProblem()->getSolution(0);
            *phaseProb_.getProblem()->getSolution(0) << function_(noCell(), X());
            //refinement_->refine(10, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
            //                                  valueOf(phaseProb_.getProblem()->getSolution(0))));
            auto &ownPhase = *phaseProb_.getProblem()->getSolution(0);
            std::cerr << ownPhase.l1norm() << "\n" << "*******";
            //double x = ((rank_ + 1) / phaseProb_.numInitialCells_) - 1;
            //int y = int(log2(x));
            //int parent = rank_ - phaseProb_.numInitialCells_ * (1 << y);

            //experimental
             
            for(int i = 0; i < phaseProb_.waitingCells_.size(); i++){
              if(rank_ == phaseProb_.waitingCells_[i]){
                daughter_index = i;
              }
            }
            //auto indexItr = std::find(phaseProb_.waitingCells_.begin(), phaseProb_.waitingCells_.end(), rank_); //finds index of daughter cell in waiting cells
            //int daughter_index = indexItr - phaseProb_.waitingCells_.begin();
            //auto daughter_index = std::distance(phaseProb_.waitingCells_.begin(), indexItr);
            int mother_rank = phaseProb_.aliveCells_[daughter_index];
            std::cerr << "\n my mother is " << mother_rank << "and I am " << rank_;
            char from_r[10];
            sprintf(from_r, "%03d", mother_rank);//phaseProb_.parent_);
            //sprintf(from_r, "%03d", phaseProb_.aliveCells_[daughter_index]);
            char to_r[10];
            sprintf(to_r, "%03d", rank_);

            Parameters::get("phase->growth_factor", phaseProb_.growth_factor_);

            //AMDiS::io::readFile(str1+str2+str3, ownPhase);
            std::string arh_filename{arh_file_prefix_ + std::string(from_r) + "_to_rank_" + std::string(to_r) + ".arh"};
            std::cout << arh_filename + "\n";
            AMDiS::io::readFile(arh_filename, ownPhase);
            
            
            //phaseProb_.calcSignedDist(adaptInfo);
            //refinement_->refine(50, function_(indicator2(getMesh()->getName() + "->phase", 3 * eps_),
            //                           valueOf(*phaseProb_.tmp_lagrange1_)));
            //refinement_->refine(300, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
            //                                  valueOf(getProblem()->getSolution(0))));
          }
        }
      }

      if(status>0){
        //refinement_->refine(100, function_(indicator3(getMesh()->getName() + "->phase"/*, 7 * eps_*/),
        //                               valueOf(*phaseProb_.getProblem()->getSolution(0))));
        
        
        phaseProb_.calcSignedDist(adaptInfo, true);
        refinement_->refine(30, function_(indicator4(getMesh()->getName() + "->phase", 7 * eps_),
                                      valueOf(*phaseProb_.tmp_lagrange1_)));
        for (int _rep = 0; _rep < 5; ++_rep)
        {
          phaseProb_.calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(3, function_(indicator4(getMesh()->getName() + "->phase", 7 * eps_),
                                          valueOf(*phaseProb_.tmp_lagrange1_)));
        }
        // calculate phase-field from signed-dist function
        *getProblem()->getSolution(0) << tanh(valueOf(*phaseProb_.tmp_lagrange1_) / (-std::sqrt(2.0) * eps_));

        ///AMDiS::Constant One(1.0);
        ///ElementFunctionAnalytic<double> oneFct(&One);
        ///ElementFunctionDOFVec<double> elLevelFct(phaseProb_.tmp_lagrange1_.get());
        ///compositeFEM::ElementLevelSet elLevelSet("dist", &elLevelFct, getProblem()->getMesh());
        phaseProb_.pos_ = integrate(X() * ((valueOf(getProblem()->getSolution(0)) + 1.0) / 2.0)) / integrate((valueOf(getProblem()->getSolution(0)) + 1.0) / 2.0);

        // phaseProb_.current_volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
        phaseProb_.current_volume_ = integrate(0.5*valueOf(getProblem()->getSolution(0))+0.5);
        
        phaseProb_.radius_ = std::sqrt(phaseProb_.current_volume_ / M_PI);
        
      }

      MPI_Barrier(comm());
      
      //rescaling of Ca_
      ///(phaseProb_.current_volume_ < volume_threshold_*0.5) ? phaseProb_.rescaledCa_ = 0.1 : phaseProb_.rescaledCa_=phaseProb_.Ca_;

      t.reset();
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
      std::vector<std::pair<std::size_t,double>> contributions;

      communicator_->gather_and_deliver(phi_j, contributions, [&phi_i, &interaction, b1_w=this->B1_wr(), w1_b=this->wr1_B()](int macroIndex, DOFVector<double> const& phi_j)
        {
          return transfer_mm(phi_i, phi_j, interaction, macroIndex, [b1_w,w1_b](double p_i, double p_j)
          {
            return b1_w(p_i,p_j) + w1_b(p_i,p_j);
          });
        });         
      MSG("time (gather-scatter) = %e\n", t.elapsed());
      MPI_Barrier(comm());

      // ----- Neighbour and Interaction output from here ------------
      std::string directory_ = ".";
      std::string postfix_ = "";
      Parameters::get("output", directory_);
      Parameters::get("postfix", postfix_);
      postfix_ += "_p" + std::to_string(comm().rank());
      
      
      bool writeAliveAndWaitingCells = true;
      if (writeAliveAndWaitingCells){
        std::string alive_filename = directory_ + "/alive" + postfix_ + ".dat";
        std::string wait_filename = directory_ + "/wait" + postfix_ + ".dat";
        std::ofstream aliveFile;
        std::ofstream waitFile;
        aliveFile.open(alive_filename,std::ios::app);
        waitFile.open(wait_filename,std::ios::app);
        for (int r = 0; r<phaseProb_.aliveCells_.size(); ++r){
          aliveFile << phaseProb_.aliveCells_[r] << " ";
          waitFile << phaseProb_.waitingCells_[r] << " ";
        }
        aliveFile <<  "\n";
        waitFile <<  "\n";
        aliveFile.close();  
        waitFile.close(); 
      }

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

      if(contact_inhi_growth == 1){
        double total_contributions;
        phaseProb_.neighbours = 0;
        for (std::size_t i = 0; i < contributions.size(); i++){
          total_contributions += contributions[i].second;
          if(abs(contributions[i].second) > 0.0){//non zero contributions means that a cell has a neighbour
            phaseProb_.neighbours += 1;
          }
        }
        total_contributions = total_contributions/In_;        
        phaseProb_.total_interactions = total_contributions;
        phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_*(std::max(0.0, std::min(1.0 - ((total_contributions > 0) - (total_contributions < 0))*std::pow((total_contributions/inhibiting_strength),2),1.0)));
        //phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_*(std::max(0.0, std::min(1.0 - (total_contributions/std::abs(total_contributions))*std::pow((total_contributions/inhibiting_strength),2),1.0)));
        //if(phaseProb_.current_volume_ < volume_threshold_*0.5){   //if the cell is small, there is no contact inhibition
        //  phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_; 
        //}
      }     

      /*if(status > 0 ){
        phaseProb_.gamma_ = 0.0; //disabling CH if cutting just s
        if(rank_==0){
        AMDiS::io::VtkWriter::writeFile(*phaseProb_.getProblem()->getSolution(0), "/scratch/ws/1/haja565a-workspace1/dat4.vtu");}
        //interaction_->set(0.0);
      } else {
        Parameters::get("phase->gamma", phaseProb_.gamma_);
      }*/

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
    std::vector<std::vector<double>> volumeOut_;
    std::vector<std::vector<double>> volumeIn_;

  };
} } // end namespaces
