#pragma once

#include <fstream>
#include <limits>
#include <cstdio>
#include <random>

#include <mpi14/Collective.hpp>

#include "InitialFunctions.h"
#include "CutPhase.h"
#include "Problem.h"
#include "Refinement.h"
#include "VerletList2.h"
#include "MPIBaseProblem.h"
#include "include/ContourWriter.h"
#include "MacroCommunicator.h"

#include "HL_SignedDistTraverse.h"
#include "compositeFEM/CFE_Integration.h"
#include "extensions/RefinementExpression.h"

#include<typeinfo>

namespace AMDiS
{
namespace base_problems
{

class PhaseFieldProblem
    : public MPIBaseProblem<Problem>
{
  using ProblemType = Problem;
  using Super = MPIBaseProblem<ProblemType>;

  friend class PhaseFieldGlobal;

public:
  PhaseFieldProblem(mpi14::Communicator comm)
      : Super("phase", comm)
  {
    Parameters::get("phase->gamma", gamma_);
    Parameters::get("phase->epsilon", eps_);
    Parameters::get("phase->In", In_);
    Parameters::get("phase->Ca", Ca_);
    Parameters::get("global->v0", v0_);
    Parameters::get("global->D", angleDiffusivity_);
    Parameters::get("phase->allen_growth", allen_growth_);
  	Parameters::get("phase->growth_factor", growth_factor_);
    Parameters::get("phase->division_volume", volume_threshold_);
    actual_growth_factor_ = growth_factor_; //changes based on interaction possible
    Parameters::get("phase->growth_type", growth_type_);
    Parameters::get("phase->Gr", Gr_);
    Parameters::get("phase->gamma2", gamma2_);
    Parameters::get("phase->So", So_);
    Parameters::get("phase->confined", confined_);
    Parameters::get("phase->Con", Con_);
    Parameters::get("phase->con_inhibition_limit", con_inhibiting_limit_);
    Parameters::get("phase->confinement_size", confinement_size_);
    Parameters::get("phase->confinement_shape", confinementShape_);
    Parameters::get("phase->switch_size", switchsize);
  }

  virtual void initData() override
  {
    // INITIALIZE ALL VECTORS
    // tmp lagrange for approximate signed distance functions...
    tmp_lagrange2_.reset(new DOFVector<double>(getFeSpace(0), "tmp_lagrange2"));
    tmp_lagrange1_.reset(new DOFVector<double>(getFeSpace(2), "tmp_lagrange1"));
    // advection is computed using a brownian motion
    advection_[0]  = new DOFVector<double>(getFeSpace(0), "velocity_x");
    advection_[1]  = new DOFVector<double>(getFeSpace(0), "velocity_y");
    // elongation tensor, used in postprocessing to evaluate the shape of cells
    nematic_tensor_[0]  = 1.0;
    nematic_tensor_[1]  = 0.0;
    major_axis_angle_ = 90.0;



    // random angle in the beginning
    int N = 100;
    Parameters::get("number of cells", N);
    N = std::max(N, 100);
    std::srand(random_seed);
    std::vector<double> rand_angles(N);
    std::generate(rand_angles.begin(), rand_angles.end(), [](){ return std::rand()/double(RAND_MAX)*2.0*M_PI; });
    theta_ = rand_angles[rank_];


    reinit_.reset(new HL_SignedDistTraverse("reinit", 2));
    refinement_.reset(new extensions::RefinementExpression(getProblem()->getMesh()));
    verletList_.reset(new VerletList(getFeSpace(0)->getMesh(), comm_));
	

    // OUTPUT DIRECTORIES
    directory_ = ".";
    postfix_ = "";
    setup_= "";
    Parameters::get("output", directory_);
    Parameters::get("postfix", postfix_);
    Parameters::get("main->setup", setup_);
    postfix_ += "_p" + std::to_string(rank_);

    std::remove((directory_ + "/data" + postfix_ + ".csv").c_str());
    std::remove((directory_ + "/interaction" + postfix_ + ".dat").c_str());
    std::remove((directory_ + "/velocities" + postfix_ + ".dat").c_str());
    std::remove((directory_ + "/neighbours" + postfix_ + ".dat").c_str());

    std::string new_filename1 = directory_ + "/signedDist" + postfix_ + "_";
    Parameters::set("signedDist->output->filename", new_filename1);

    filename_positions_ = directory_ + "/positions" + postfix_ + ".csv";
    writer_signedDist_.reset(new ContourWriter("signedDist->output", tmp_lagrange1_.get()));

    writeData_ = false;
    Parameters::get("main->write positions", writeData_);
    Parameters::get(getMesh()->getName() + "->dimension",domainDimension_);
    if (writeData_) {
        std::ofstream out(filename_positions_);
        out << "time,rank,x0,x1,r,S0,S1,v0,v1,angle, total_interaction, neighbours, confine_interaction, growth_rate\n";
    }

    if (confined_==1){
      //domainDimension_ = rescaleConfinement(domainDimension_); 

      confinementPF_.reset(new DOFVector<double>(getFeSpace(0), "confinement phase field"));
      confinementForce_.reset(new DOFVector<double>(getFeSpace(0), "confinement force"));

      if((confinementShape_ == 1)){
        *confinementPF_ << function_(invertedRectangle(0.5*domainDimension_, eps_, domainDimension_[0], confinement_size_*domainDimension_[0], confinement_size_*domainDimension_[1], false), X());
      }else if((confinementShape_ == 2)){
        *confinementPF_ << function_(invertedRotatedEllipse(0.5*domainDimension_, eps_, domainDimension_[0], 0.0, 0.5*confinement_size_*domainDimension_[0], 0.5*confinement_size_*domainDimension_[1], false), X());
      }else if((confinementShape_ == 3)){
        *confinementPF_ << function_(annularRing(0.5*domainDimension_, eps_, domainDimension_[0], 0.0, confinement_size_*domainDimension_[0], 0.10*domainDimension_[1], false), X());
      }
      else
        TEST_EXIT(false)("Undefined confinement shape.");  
    }

  }

  virtual void solveInitialProblem(AdaptInfo *adaptInfo) override
  {
    FUNCNAME("PhaseFieldProblem::solveInitialProblem()");
    MSG("start (PhaseFieldProblem::solveInitialProblem) = %20.16e\n", mpi14::now());
    bool gaussianSize = 0;
    Parameters::get(name_ + "->initial->gaussian", gaussianSize);

    if (setup_.compare("singular") == 0)
      singularInitialCell(adaptInfo);
    else if (setup_.compare("collective") == 0){
      if (gaussianSize)
        gaussianRectangles(adaptInfo);
      else
        hexagonalRectangleInitial(adaptInfo);
    } 
    else if(setup_.compare("binary_Lea")==0) {
      binaryLeaInitial(adaptInfo);
    }
    else if(setup_.compare("three_collision")==0){
      threeCollisionInitial(adaptInfo);
    }
    else if(setup_.compare("two_division")==0){
      TwoDivisionInitial(adaptInfo);
    }
    else if(setup_.compare("one_division")==0){
      OneDivisionInitial(adaptInfo);			
    }
    else if(setup_.compare("multi_cell")==0){
      MultiCellGrowDivide(adaptInfo);
    }
    else if(setup_.compare("multi_ring")==0){
      MultiCellMultiRing(adaptInfo);
    }
    else
      TEST_EXIT(false)("Undefined setup.");  

    calcSignedDist(adaptInfo, true);

    // calculate initial phase-field from signed-dist function
    *getProblem()->getSolution(0) << tanh(valueOf(*tmp_lagrange1_) / (-std::sqrt(2.0) * eps_));

    pos_= integrate(X() * ((valueOf(getProblem()->getSolution(0)) + 1.0) / 2.0)) / integrate((valueOf(getProblem()->getSolution(0)) + 1.0) / 2.0);

    initial_volume_ = integrate((valueOf(getProblem()->getSolution(0)) + 1.0) / 2.0);
    
    radius_ = std::sqrt(initial_volume_/M_PI);
    current_volume_ = initial_volume_;

    // calculate grid-width
    h_ = mesh_width(getMesh());

    MSG("num DOFs[1]: %d\n", tmp_lagrange1_->getFeSpace()->getAdmin()->getUsedDofs());
    MSG("num DOFs[2]: %d\n", tmp_lagrange2_->getFeSpace()->getAdmin()->getUsedDofs());
    MSG("mesh width: %f\n", h_);

    refinement_->setRefineOperation(adaptInfo, getProblem());
  }

  void MultiCellGrowDivide(AdaptInfo *adaptInfo){
    std::cout << "initialise some cells and rankTrees" << std::endl;
    
    Parameters::get(name_ + "initial->resize", resize);
    resize = 0.25;
    //create rank tree
    Parameters::get(name_ + "->initial->numCells", numInitialCells_);
    int treesize = size_ / numInitialCells_ + 1;
    rankTree_.resize(treesize);
    int factor = 0;
    alive_ = 1;

    if (rank_ >= numInitialCells_)
    {
      alive_ = 0; 
      double x = rank_ / numInitialCells_;
      factor = int(log2(x)) + 1;
    }
    rankTree_[0] = rank_;
    for (int j = 1; j < treesize; ++j)
    {
      rankTree_[j] = rank_ + ((1 << (j + factor - 1)) * numInitialCells_); //x<<n == x*(2^n)
    }
    nextChildIndex_ = 1;
    parent_ = rank_; //parent of initial cells are themselves(dummy). 

    if (rank_ >= numInitialCells_)
    {
      double x = ((rank_ + 1) / numInitialCells_) - 1;
      int y = int(log2(x));
      parent_ = rank_ - numInitialCells_ * (1 << y);
    }

    std::cout << rankTree_ << "\n";
    
    //experimental ranking of cells
    aliveCells_.resize(numInitialCells_);
    waitingCells_.resize(numInitialCells_);
    latestWaitingCell_ = numInitialCells_*2-1;
    std::iota(aliveCells_.begin(), aliveCells_.end(), 0);
    std::iota(waitingCells_.begin(), waitingCells_.end(), numInitialCells_);
    //std::cerr << aliveCells_ << "\n";
    //std::cerr << waitingCells_ << "\n";
    //std::cerr << latestWaitingCell_ << "\n";

    //Improve this later to suit more cells
    WorldVector<double> position;

    if(numInitialCells_==1){
      if (rank_ == 0)
      {
        alive_ = 1;
        position[0] = 0.5 * domainDimension_[0];
        position[1] = 0.5 * domainDimension_[1];
      
        double a = 0.25 * resize * domainDimension_[0];
        double b = 0.25 * resize * domainDimension_[1];

        *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

        for (int _rep = 0; _rep < 5; ++_rep)
        {
          calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                           valueOf(getProblem()->getSolution(0))));

          *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        }
      }
      else if (rank_ >= 1)
      {
        position[0] = 0.0;
        position[1] = 0.0;
        double a = 0.01 * resize * domainDimension_[0];
        double b = 0.01 * resize * domainDimension_[1];
        //*getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        //growth_factor_ = 0.0;
        *getProblem()->getSolution(0) << function_(noCell(), X());
      }
    }
    if(numInitialCells_ > 1){
      if (rank_ == 0)
      {
        alive_ = 1;
        position[0] = 0.5 * domainDimension_[0];
        position[1] = 0.5 * domainDimension_[1];

        //set size for collision
        int d = 1.0;
        if(switchsize==2){ d = 0.0;}
      
        double a = 0.25 * resize * domainDimension_[0]*d;
        double b = 0.25 * resize * domainDimension_[1]*d;

        *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

        for (int _rep = 0; _rep < 5; ++_rep)
        {
          calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                           valueOf(getProblem()->getSolution(0))));

          *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        }
      }
      
      double angle_difference_ = 360.0 / (numInitialCells_-1);
      if(rank_<numInitialCells_ && rank_ >0){
        double orientation_ = (rank_-1) * angle_difference_;
        alive_ = 1;
        position[0] = domainDimension_[0]*(0.5 + 0.15*std::cos(orientation_ * M_PI / 180.0));
        position[1] = domainDimension_[1]*(0.5 + 0.15*std::sin(orientation_ * M_PI / 180.0));
        
        //growth_factor_ = 0.0;
        
        //size for cell collision
        int c = 1.0;
        if(switchsize==1){ c = 0.0;}

        double a = 0.25 * resize * domainDimension_[0]*c;
        double b = 0.25 * resize * domainDimension_[1]*c;

        if (std::max(a, b) * numInitialCells_ * 2 > 6 * 0.25 * domainDimension_[0])
        {
          std::cerr << "Too many cells " << std::endl;
          std::abort();
        }

        *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

        for (int _rep = 0; _rep < 5; ++_rep)
        {
          calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                           valueOf(getProblem()->getSolution(0))));

          *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        }
      }
      else if(rank_ >= numInitialCells_){
        position[0] = 0.0;
        position[1] = 0.0;
        alive_ = 0;
        double a = 0.01 * resize * domainDimension_[0];
        double b = 0.01 * resize * domainDimension_[1];
        //*getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        //growth_factor_ = 0.0;
        *getProblem()->getSolution(0) << function_(noCell(), X());
      }
    }
  }
  void MultiCellMultiRing(AdaptInfo *adaptInfo){
    std::cout << "initialise some cells and rankTrees" << std::endl;
    
    Parameters::get(name_ + "initial->resize", resize);
    //create rank tree
    Parameters::get(name_ + "->initial->numCells", numInitialCells_);
    int treesize = size_ / numInitialCells_ + 1;
    rankTree_.resize(treesize);
    int factor = 0;
    alive_ = 1;
    if (rank_ >= numInitialCells_)
    {
      alive_ = 0; 
      double x = rank_ / numInitialCells_;
      factor = int(log2(x)) + 1;
    }
    rankTree_[0] = rank_;
    for (int j = 1; j < treesize; ++j)
    {
      rankTree_[j] = rank_ + ((1 << (j + factor - 1)) * numInitialCells_); //x<<n == x*(2^n)
    }
    nextChildIndex_ = 1;
    parent_ = rank_; //parent of initial cells are themselves(dummy). 

    if (rank_ >= numInitialCells_)
    {
      double x = ((rank_ + 1) / numInitialCells_) - 1;
      int y = int(log2(x));
      parent_ = rank_ - numInitialCells_ * (1 << y);
    }

    std::cout << rankTree_ << "\n";
    
    //experimental ranking of cells
    aliveCells_.resize(numInitialCells_);
    waitingCells_.resize(numInitialCells_);
    latestWaitingCell_ = numInitialCells_*2-1;
    std::iota(aliveCells_.begin(), aliveCells_.end(), 0);
    std::iota(waitingCells_.begin(), waitingCells_.end(), numInitialCells_);
    //std::cerr << aliveCells_ << "\n";
    //std::cerr << waitingCells_ << "\n";
    //std::cerr << latestWaitingCell_ << "\n";

    //Improve this later to suit more cells
    WorldVector<double> position;

    
      if (rank_ == 0)
      {
        alive_ = 1;
        position[0] = 0.5 * domainDimension_[0];
        position[1] = 0.5 * domainDimension_[1];
      
        double a = 0.25 * resize * domainDimension_[0];
        double b = 0.25 * resize * domainDimension_[1];

        *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

        for (int _rep = 0; _rep < 5; ++_rep)
        {
          calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                           valueOf(getProblem()->getSolution(0))));

          *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        }
      }
      
      double angle_difference_ = 360.0 / (numInitialCells_-1);
      if(rank_<numInitialCells_ && rank_ >0){
        double orientation_ = (rank_-1) * angle_difference_;
        alive_ = 1;

        if(rank_%5==2){
          position[0] = domainDimension_[0]*(0.5 + 0.15*std::cos(orientation_ * M_PI / 180.0));
          position[1] = domainDimension_[1]*(0.5 + 0.15*std::sin(orientation_ * M_PI / 180.0));
        }else if(rank_%5==0 || rank_%5==3){
          position[0] = domainDimension_[0]*(0.5 + 0.42*std::cos(orientation_ * M_PI / 180.0));
          position[1] = domainDimension_[1]*(0.5 + 0.42*std::sin(orientation_ * M_PI / 180.0));
        }else{
          position[0] = domainDimension_[0]*(0.5 + 0.27*std::cos(orientation_ * M_PI / 180.0));
          position[1] = domainDimension_[1]*(0.5 + 0.27*std::sin(orientation_ * M_PI / 180.0));
        }
        //growth_factor_ = 0.0;

        double a = 0.23 * resize * domainDimension_[0];
        double b = 0.23 * resize * domainDimension_[1];

        /*if (std::max(a, b) * numInitialCells_ * 2 > 6 * 0.25 * domainDimension_[0])
        {
          std::cerr << "Too many cells " << std::endl;
          std::abort();
        }*/

        *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

        for (int _rep = 0; _rep < 5; ++_rep)
        {
          calcSignedDist(adaptInfo);

          // refine/coarsen local mesh along interface
          refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                           valueOf(getProblem()->getSolution(0))));

          *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        }
      }
      else if(rank_ >= numInitialCells_){
        position[0] = 0.0;
        position[1] = 0.0;
        alive_ = 0;
        double a = 0.01 * resize * domainDimension_[0];
        double b = 0.01 * resize * domainDimension_[1];
        //*getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
        //growth_factor_ = 0.0;
        *getProblem()->getSolution(0) << function_(noCell(), X());
      }
    
  }
  void singularInitialCell(AdaptInfo *adaptInfo){
    /*
    * Create single cell for activity treshold experiment
    * Initial center is the center of the domain
    * Initial diameter is dimension*resize
    */
    std::cout << "enter singularInitialCell" << std::endl;
    double resize = 1.0;
    Parameters::get(name_ + "->initial->resize", resize);

    WorldVector<double> position;
    for (int i = 0; i < dow_; i++)
      position[i] = 0.5 * domainDimension_[i];

    double a = 0.5 * resize * domainDimension_[0];
    double b = 0.5 * resize * domainDimension_[1];
    std::cout << "a" << a << std::endl;
    std::cout << "b" << b << std::endl;
    std::cout << "positoio" << position[0] << position[1]<<std::endl;
    *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

    for (int _rep = 0; _rep < 5; ++_rep)
    {
      calcSignedDist(adaptInfo);

      // refine/coarsen local mesh along interface
      refinement_->refine(5, function_(indicator2(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(*tmp_lagrange1_)));

      *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

    }
  }
  void binaryLeaInitial(AdaptInfo *adaptInfo){
    /*
    Creates two cells (aligned on x-axis, but not aligned on y-axis). 
    Size is determined by dimension and rescaled.
    */
   std::cout << "enter binaryLeaInitial" << std::endl;
   double resize=1.0;
   Parameters::get(name_+"->initial->resize",resize);
   WorldVector<double> position;
   position[0]=0.5*domainDimension_[0];
   position[1]=(0.25+rank_ *0.5)*domainDimension_[1];
   double a= 0.5*resize*domainDimension_[0];
   double b=0.5*resize*domainDimension_[1];
   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "position: (" << position[0] << " , " << position[1] << ")" << std::endl;
   *getProblem()->getSolution(0) << function_(rectangle(position,eps_,domainDimension_[0],a,b),X());
   for (int _rep = 0; _rep < 5; ++_rep)
    {
      calcSignedDist(adaptInfo);

      // refine/coarsen local mesh along interface
      refinement_->refine(5, function_(indicator2(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(*tmp_lagrange1_)));

      *getProblem()->getSolution(0) << function_(rectangle(position, eps_, domainDimension_[0], a, b), X()); 
    }
  }
  void threeCollisionInitial(AdaptInfo *adaptInfo){
    /* Creates 3 cells, start with circle
    */
   std::cout << "enter threeCollisionInitial" << std::endl;
   double resize=1.0;
   Parameters::get(name_+"->initial->resize",resize);
   WorldVector<double> position;
   if (rank_==0){
    position[0]=0.5*domainDimension_[0];
    position[1]=0.75*domainDimension_[1];
   }
   else if (rank_==1){
    position[0]=0.25*domainDimension_[0];
    position[1]=(0.75-std::sqrt(3)/4.0)*domainDimension_[1];
   }
   else{
    position[0]=0.75*domainDimension_[0];
    position[1]=(0.75-std::sqrt(3)/4.0)*domainDimension_[1];
   }
   double a= 0.25*resize*domainDimension_[0];
   double b=0.25*resize*domainDimension_[1];
    *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

    for (int _rep = 0; _rep < 5; ++_rep)
    {
      calcSignedDist(adaptInfo);

      // refine/coarsen local mesh along interface
      refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(getProblem()->getSolution(0))));

      *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

    }
  }
  
  void TwoDivisionInitial(AdaptInfo *adaptInfo){
    /* Creates 2 cells and also 2 empty phases to accomodate newly divided phases, start with circle
    */
   std::cout << "enter twoDivisionInitial" << std::endl;
   double resize=1.0;
   Parameters::get(name_+"->initial->resize",resize);
   WorldVector<double> position;
   if (rank_==0){
    position[0]=0.5*domainDimension_[0];
    position[1]=0.75*domainDimension_[1];
   }
   else if (rank_==1){
    position[0]=0.5*domainDimension_[0];
    position[1]=0.25*domainDimension_[1];
   }
   if (rank_ <= 1){
	   double a= 0.25*resize*domainDimension_[0];
	   double b=0.25*resize*domainDimension_[1];
       *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

       for (int _rep = 0; _rep < 5; ++_rep)
       {
			calcSignedDist(adaptInfo);

			// refine/coarsen local mesh along interface
			refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(getProblem()->getSolution(0))));

			*getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
		}
   }
   else if(rank_ > 1){
	   *getProblem()->getSolution(0) << function_(noCell(), X());
   }
  }
  
  void OneDivisionInitial(AdaptInfo *adaptInfo){
    /* Creates 1 cells and also 3 empty phases to accomodate newly divided phases, start with circle
    */
   std::cout << "enter oneDivisionInitial" << std::endl;
   double resize=1.0;
   Parameters::get(name_+"->initial->resize",resize);
   WorldVector<double> position;
   if (rank_==0){
    position[0]=0.5*domainDimension_[0];
    position[1]=0.5*domainDimension_[1];
   }
   if (rank_ < 1){
	   double a= 0.25*resize*domainDimension_[0];
	   double b=0.25*resize*domainDimension_[1];
       *getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());

       for (int _rep = 0; _rep < 5; ++_rep)
       {
			calcSignedDist(adaptInfo);

			// refine/coarsen local mesh along interface
			refinement_->refine(5, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(getProblem()->getSolution(0))));

			*getProblem()->getSolution(0) << function_(rotatedEllipse(position, eps_, domainDimension_[0], 0.0, a, b, false), X());
		}
   }
   else if(rank_ >= 1){
	   *getProblem()->getSolution(0) << function_(noCell(), X());
   }
  }
  
  void hexagonalRectangleInitial(AdaptInfo *adaptInfo){
    /*
    Create dense rectangles in hexagonal alignment. Size is determined by 
    dimension and rescaled.
    */
    double resize = 1.0;
    Parameters::get(name_ + "->initial->resize", resize);

    int num_cells = 1;
    Parameters::get("number of cells", num_cells);

    int Nx = std::floor(std::sqrt(num_cells));

    int column = rank_ % Nx; //column
    int row = rank_ / Nx;     //row

    double sizeX = domainDimension_[0] / Nx; //size in dim 0
    double sizeY = domainDimension_[1] / Nx; //size in dim 1

    WorldVector<double> position;
    position[0] = (row % 2 == 0 ? 0.5 * sizeX + column * sizeX : column * sizeX);

    position[1] = 0.5 * sizeY + row * sizeY;

    sizeX *= resize;
    sizeY *= resize;

    *getProblem()->getSolution(0) << function_(rectangle(position, eps_, domainDimension_[0], sizeX, sizeY), X());
   
    for (int _rep = 0; _rep < 5; ++_rep)
    {
      calcSignedDist(adaptInfo);

      // refine/coarsen local mesh along interface
      refinement_->refine(5, function_(indicator2(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(*tmp_lagrange1_)));

      *getProblem()->getSolution(0) << function_(rectangle(position, eps_, domainDimension_[0], sizeX, sizeY), X()); 
    }
  }

  void gaussianRectangles(AdaptInfo *adaptInfo){
    /*
    * Constant row height but inside the rows a normal distribution of cell sizes
    * is computed. The last cell computes this and then adjusts its own size.
    */
    double resize = 1.0;
    Parameters::get(name_ + "->initial->resize", resize);

    int num_cells = 1;
    Parameters::get("number of cells", num_cells);

    int Nx = std::floor(std::sqrt(num_cells));

    int column = rank_ % Nx; //column
    int row = rank_ / Nx;     //row    

    double sizeY = domainDimension_[1] / Nx; //size in dim 1
    // ... until here the same as for equidistant rectangles

    // only declared but later set, possibly from other rank
    double sizeX = 0.0; 

    // only for y direction the position is clear
    WorldVector<double> position;
    position[1] = 0.5 * sizeY + row * sizeY;

    // Size in x direction is no longer deterministic
    double meanSize = domainDimension_[0] / Nx;   

    if (column == Nx - 1){
      // Last rank in a row ... do the computation

      // We need to loop until we have a 'good' distribution, i.e.
      // no negative values 
      //bool check = true;
      std::vector<double> sizes;

      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::normal_distribution<double> distribution(meanSize,meanSize/4.0);
        
      for (int i = 0; i < Nx; i++)
        sizes.push_back(distribution(generator));
      
      // Now we normalize
      double size_sum = std::accumulate(sizes.begin(), sizes.end(), 0.0);
      double norm_factor = domainDimension_[0]/size_sum;

      for (int i = 0; i < Nx; i++)
        sizes[i] *= norm_factor;

      // Here we assign the value for later ...
      sizeX = sizes[Nx-1];

      // Now compute the x-positions for all ranks
      std::vector<double> xPositions;
      for (int i = 0; i < Nx; i++){
        xPositions.push_back(row%2==0?0.5*sizes[0]:0.0);
        
        if (i>0){
          xPositions[i] += 0.5*sizes[0];
          xPositions[i] += std::accumulate(sizes.begin()+1, sizes.begin()+i, 0.0);
          xPositions[i] += 0.5 * sizes[i];
        }
      }  
      std::cout << sizes << std::endl;
      std::cout << xPositions << std::endl;
      // Now send stuff ...
      for (int i = 0; i < Nx-1; i++){
        comm_.send(sizes[i], row*Nx+i);
        comm_.send(xPositions[i], row*Nx+i);
      }

      // Here we assign the value for later..
      position[0] = xPositions[Nx-1];
    }
    else {
      comm_.recv(sizeX, (row+1)*Nx-1);
      comm_.recv(position[0], (row+1)*Nx-1);
    }

    sizeX *= resize;
    sizeY *= resize;
    *getProblem()->getSolution(0) << function_(rectangle(position, eps_, domainDimension_[0], sizeX, sizeY), X());
   
    for (int _rep = 0; _rep < 5; ++_rep)
    {
      calcSignedDist(adaptInfo);

      // refine/coarsen local mesh along interface
      refinement_->refine(5, function_(indicator2(getMesh()->getName() + "->phase", 3 * eps_),
                                       valueOf(*tmp_lagrange1_)));

      *getProblem()->getSolution(0) << function_(rectangle(position, eps_, domainDimension_[0], sizeX, sizeY), X()); 
    }
    
  }

  void calcSignedDist(AdaptInfo *adaptInfo, bool updateVerletList = false)
  {
    FUNCNAME("PhaseFieldProblem::calcSignedDist()");
    MSG("start (PhaseFieldProblem::calcSignedDist) = %20.16e\n", mpi14::now());

    // reinit phase-field
    Timer t;
    tmp_lagrange1_->interpol(getProblem()->getSolution(0), -1.0); // copy to lagrange 1 feSpace
    reinit_->calcSignedDistFct(adaptInfo, tmp_lagrange1_.get());  // redistancing
    MSG("time (reinit) = %e\n", t.elapsed());

    if (updateVerletList)
    {
      t.reset();
#ifdef DEBUG_PERFORMANCE
      MPI_Barrier(comm_);
#endif
      MSG("time (barrier) = %e\n", t.elapsed());

      writer_signedDist_->writeFiles(adaptInfo, false);
      // io::VtkWriter::writeFile(tmp_lagrange1_.get(), directory_ + "/signedDist_" + std::to_string(adaptInfo->getTimestepNumber()) +  ".vtu");

      t.reset();
      double dist = 6 * eps_;
      Parameters::get("phase->interaction", dist);
      verletList_->compute(adaptInfo, *tmp_lagrange1_, dist);
	    //verletList_->computeNoOverlap(adaptInfo, *tmp_lagrange1_, dist);
      MSG("time (updateVerletList) = %e\n", t.elapsed());
    }
  }

  // B'(phi_i) * w(phi_j)
  auto B1_w(DOFVector<double> *phi_i_vec, DOFVector<double> *phi_j_vec)
  {
    return function_([eps = eps_](double phi_i, double phi_j) {

      phi_i = std::max(-1.0, std::min(1.0, phi_i));
      phi_j = std::max(-1.0, std::min(1.0, phi_j));
      double psi_i = (phi_i-1.0)/2;
      double psi_j = (phi_j-1.0)/2;
      double a = 1.0; //coefficient of x^4
      double b = -2.0; //coefficient of x^2
      double B_der = 0.5;
      double w = 1 + b*sqr(psi_j) + a*sqr(sqr(psi_j));
      return phi_i > -0.99 && phi_j > -0.99
              ? B_der*w
              : 0.0;
    }, valueOf(phi_i_vec), valueOf(phi_j_vec));
  }

  virtual void initTimestep(AdaptInfo *adaptInfo) override
  {
    FUNCNAME("PhaseFieldProblem::initTimestep()");
    MSG("start (PhaseFieldProblem::initTimestep) = %20.16e\n", mpi14::now());


    Super::initTimestep(adaptInfo);
    if (setup_.compare("binary_Lea")==0 && adaptInfo->getTime()<=8.0){
    *advection_[0] << v0_*0.0;
    *advection_[1] << v0_*(-1.0+rank_*2.0);
    }
    else if(setup_.compare("three_collision")==0 && adaptInfo->getTime()<=8.0){
      if(rank_==0){
        *advection_[0]<<v0_*0.0;
        *advection_[1]<<v0_*(1.0);
      }
      else if(rank_==1){
        *advection_[0] << -v0_*std::sqrt(3)*0.5;
        *advection_[1] << -v0_*0.5;
      }
      else{
        *advection_[0] << v0_*std::sqrt(3)*0.5;
        *advection_[1] << -v0_*0.5;
      }
    }
    else{
    // Random advection term, Wiener process w.r.t old direction of advection
    random_seed += 1;
    std::srand(random_seed); //remove this
    std::random_device rd;
    std::mt19937 generator(random_seed);//std::mt19937 generator(rd()); 

    std::normal_distribution<double> distribution(0,1);    
    
    theta_ += angleDiffusivity_ * distribution(generator);
    ct = std::cos(theta_);
    st = std::sin(theta_);
    *advection_[0] << v0_ * std::cos(theta_) * (0.5*valueOf(getProblem()->getSolution(0)) + 0.5);
    *advection_[1] << v0_ * std::sin(theta_) * (0.5*valueOf(getProblem()->getSolution(0)) + 0.5);
    }
    // Compute elongation
    nematic_tensor_[0] = 0.125 * integrate(derivativeOf(getProblem()->getSolution(0),1) * derivativeOf(getProblem()->getSolution(0),1)); 
    nematic_tensor_[0] += -0.125 * integrate(derivativeOf(getProblem()->getSolution(0),0) * derivativeOf(getProblem()->getSolution(0),0));

    nematic_tensor_[1] = (-0.25) * integrate(derivativeOf(getProblem()->getSolution(0),0) * derivativeOf(getProblem()->getSolution(0),1));

    // normalize
    double norm_q = std::sqrt(2.0 * nematic_tensor_[0] * nematic_tensor_[0] + 2 * nematic_tensor_[1] * nematic_tensor_[1]);
    nematic_tensor_[0] *= (1.0 / norm_q);
    nematic_tensor_[1] *= (1.0 / norm_q);

    major_axis_angle_ = atan2(nematic_tensor_[1], nematic_tensor_[0]) * 180.0 / (2*M_PI);

    Timer t;
    // calculate phase-field from lagrange[1] signed-distance function
    ///***tmp_lagrange2_ << tanh(valueOf(*tmp_lagrange1_) / (-std::sqrt(2.0) * eps_)); // interpolate function back to lagrange 2 feSpace

    /*
    AMDiS::Constant One(1.0);
    ElementFunctionAnalytic<double> oneFct(&One);
    ElementFunctionDOFVec<double> elLevelFct(tmp_lagrange1_.get());
    compositeFEM::ElementLevelSet elLevelSet("dist", &elLevelFct, getProblem()->getMesh(0));

    volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
    double interface = compositeFEM::CFE_Integration::integrate_onZeroLs(&oneFct, &elLevelSet);
    */
    volume_ = integrate(0.5*valueOf(*getProblem()->getSolution(0))+0.5);
    if(confined_ == 1){
        double conForce = integrate(valueOf(*confinementForce_));
        confinement_inhibition_ = (std::max(0.0, 1.0 - std::pow((conForce/(Con_*con_inhibiting_limit_)),2)));
        //if(phaseProb_.current_volume_ < volume_threshold_*0.5){   //if the cell is small, there is no contact inhibition
        //  phaseProb_.actual_growth_factor_ = phaseProb_.growth_factor_; 
        //}
      }

    //***remove these comments***std::random_device rd;
    //***remove these comments***std::mt19937 generator(rd()); 
    //***remove these comments***std::normal_distribution<double> distribution(0.0,actual_growth_factor_); //no getTau here
    
    
    //***remove these comments***double del_gr_change = confinement_inhibition_*(actual_growth_factor_+ distribution(generator));// with inhibition
    //double del_gr_change = confinement_inhibition_*(growth_factor_+ distribution(generator));// no inhibition
    //***remove these comments***double del_volume = del_gr_change*volume_;
    //f_growth_ = (1.0/Gr_)*del_volume;
    //***remove these comments***f_growth_ = del_gr_change/Gr_;
    //std::cerr << f_growth_ << "\n";
    f_growth_ = 0.0;
    ///double correction = 0;// =(current_volume_ - volume_) / interface;
    
    ///*tmp_lagrange1_ << valueOf(*tmp_lagrange1_) - correction;
    //volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
    volume_ = integrate(0.5*valueOf(*getProblem()->getSolution(0))+0.5);
    //double x = 0.0;//
    double x = 0.0;// std::tanh(correction / (std::sqrt(2.0) * eps_));
    *getProblem()->getOldSolution() << valueOf(getProblem()->getSolution(0));
    ///***getProblem()->getOldSolution() << function_([x](double const &phi) {
    ///**  return std::max(-1.0, std::min(1.0, std::abs(1.0 + phi * x) > 1.e-5 ? (x + phi) / (1.0 + phi * x) : phi));
    ///**},valueOf(getProblem()->getSolution(0)));
    radius_ = std::sqrt(volume_/M_PI);
    /*
    if(growth_type_==1){
      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::normal_distribution<double> distribution(0.0,actual_growth_factor_ * (*getTau()));
        
      current_volume_ = (1.0 + confinement_inhibition_*(actual_growth_factor_ * (*getTau()))  + distribution(generator))* current_volume_;// 
      // volume measuring
      //double correction = 0.0;
      double correction;// =(current_volume_ - volume_) / interface;
      phaseStatus != 0 ? correction = 0.0: correction =(current_volume_ - volume_) / interface;
      
      *tmp_lagrange1_ << valueOf(*tmp_lagrange1_) - correction;
      volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
      MSG("volume-error = %e\n", std::abs(current_volume_ - volume_));
      
    
      //double x = 0.0;//
      double x = std::tanh(correction / (std::sqrt(2.0) * eps_));
      *getProblem()->getOldSolution() << function_([x](double const &phi) {
        return std::max(-1.0, std::min(1.0, std::abs(1.0 + phi * x) > 1.e-5 ? (x + phi) / (1.0 + phi * x) : phi));
      },
                                                  valueOf(getProblem()->getSolution(0)));
      
      radius_ = std::sqrt(volume_/M_PI);
    }else if(growth_type_==2){//this is bullshit
      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::normal_distribution<double> distribution(0.0,actual_growth_factor_); //no getTau here
      
      double gr_change = 1+actual_growth_factor_+ distribution(generator);
      double volume_new = gr_change*volume_;
      f_growth_ = (1.0/Gr_)*(1.0/volume_new)*(1.0/gr_change-1.0);
      //std::cerr << f_growth_ << "\n";

      double correction = 0;// =(current_volume_ - volume_) / interface;
      
      *tmp_lagrange1_ << valueOf(*tmp_lagrange1_) - correction;
      volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
      
      //double x = 0.0;//
      double x = std::tanh(correction / (std::sqrt(2.0) * eps_));
      *getProblem()->getOldSolution() << function_([x](double const &phi) {
        return std::max(-1.0, std::min(1.0, std::abs(1.0 + phi * x) > 1.e-5 ? (x + phi) / (1.0 + phi * x) : phi));
      },
                                                  valueOf(getProblem()->getSolution(0)));
      radius_ = std::sqrt(volume_/M_PI);
    }else if(growth_type_==3){
      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::normal_distribution<double> distribution(0.0,actual_growth_factor_); //no getTau here
      
      double del_gr_change = confinement_inhibition_*(actual_growth_factor_+ distribution(generator));//
      double del_volume = del_gr_change*volume_;
      //f_growth_ = (1.0/Gr_)*del_volume;
      f_growth_ = del_gr_change/Gr_;
      //std::cerr << f_growth_ << "\n";

      double correction = 0;// =(current_volume_ - volume_) / interface;
      
      *tmp_lagrange1_ << valueOf(*tmp_lagrange1_) - correction;
      //volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
      volume_ = integrate(0.5*valueOf(*getProblem()->getSolution(0))+0.5);
      //double x = 0.0;//
      double x = 0.0;// std::tanh(correction / (std::sqrt(2.0) * eps_));
      *getProblem()->getOldSolution() << function_([x](double const &phi) {
        return std::max(-1.0, std::min(1.0, std::abs(1.0 + phi * x) > 1.e-5 ? (x + phi) / (1.0 + phi * x) : phi));
      },
                                                  valueOf(getProblem()->getSolution(0)));
      radius_ = std::sqrt(volume_/M_PI);
    }else if(growth_type_ == 4){//hybrid growth
      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::normal_distribution<double> distribution(0.0,actual_growth_factor_);
        
      current_volume_ = (1.0 + actual_growth_factor_ * (*getTau()) + distribution(generator) * (*getTau()) ) * current_volume_;
      // volume measuring
      //double correction = 0.0;
      double correction;// =(current_volume_ - volume_) / interface;
      phaseStatus != 0 ? correction = 0.0: correction =(current_volume_ - volume_) / interface;
      
      *tmp_lagrange1_ << valueOf(*tmp_lagrange1_) - correction;
      volume_ = compositeFEM::CFE_Integration::integrate_onNegLs(&oneFct, &elLevelSet);
      MSG("volume-error = %e\n", std::abs(current_volume_ - volume_));
      
    
      //double x = 0.0;//
      double x = std::tanh(correction / (std::sqrt(2.0) * eps_));
      *getProblem()->getOldSolution() << function_([x](double const &phi) {
        return std::max(-1.0, std::min(1.0, std::abs(1.0 + phi * x) > 1.e-5 ? (x + phi) / (1.0 + phi * x) : phi));
      },
                                                  valueOf(getProblem()->getSolution(0)));
      
      radius_ = std::sqrt(volume_/M_PI);

      double del_gr_change = actual_growth_factor_+ distribution(generator);
      double del_volume = del_gr_change*volume_;
      f_growth_ = (1.0/Gr_)*del_volume;
      //std::cerr << f_growth_ << "\n";
    }*/


    bool writeVolume = false;
    Parameters::get("main->write volume", writeVolume);

    if (writeVolume)
    {
      std::ofstream out(directory_ + "/data" + postfix_ + ".csv", std::ios_base::app);
      out << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::scientific;
      out << adaptInfo->getTime() << ' ' << volume_ << '\n';
    }

    MSG("time (volume) = %e\n", t.elapsed());

    if (confined_==1){
      // Compute interaction force with confinement
      confinementForce_->set(0.0);
      *confinementForce_ << B1_w(getProblem()->getOldSolution(),confinementPF_.get());
      
      if(adaptInfo->getTimestepNumber()%50 == 0){
        std::cerr  << adaptInfo->getTimestepNumber() << ": " << confinement_interaction_*Con_ << "\n";
      }
      if(confinement_interaction_ > 0){
        if((confinementShape_ == 1)){
        *confinementPF_ << function_(invertedRectangle(0.5*domainDimension_, eps_, domainDimension_[0], confinement_size_*domainDimension_[0], confinement_size_*domainDimension_[1], false), X());
        }else if((confinementShape_ == 2)){
        *confinementPF_ << function_(invertedRotatedEllipse(0.5*domainDimension_, eps_, domainDimension_[0], 0.0, 0.5*confinement_size_*domainDimension_[0], 0.5*confinement_size_*domainDimension_[1], false), X());
        }else if((confinementShape_ == 3)){
        *confinementPF_ << function_(annularRing(0.5*domainDimension_, eps_, domainDimension_[0], 0.0, confinement_size_*domainDimension_[0], 0.10*domainDimension_[1], false), X());
        }
      }
    }
  }

  virtual void closeTimestep(AdaptInfo *adaptInfo) override
  {
    FUNCNAME("PhaseFieldProblem::closeTimestep()");
    MSG("start (PhaseFieldProblem::closeTimestep) = %20.16e\n", mpi14::now());

    calcSignedDist(adaptInfo, true);

    Timer t;
    // refine/coarsen local mesh along interface
    //refinement_->refine(1, function_(indicator2(getMesh()->getName() + "->phase", 7 * eps_),
    //                                 valueOf(*tmp_lagrange1_)));
    
    //if((current_volume_ > volume_threshold_*0.997) || div_refine_>0){ //0.49 works for volume 125.0 For higher volume, might need some change. Not yet tested
    if(div_refine_ > 0){
      boost_refinement_ = true;
      div_refine_ -= 1;
    } else{
      boost_refinement_ = false; 
    }
    
    refinement_->refine(1, function_(indicator2b(getMesh()->getName() + "->phase", 3 * eps_, boost_refinement_),
                                      valueOf(*tmp_lagrange1_)));

    ///time_counter_++;
    /*refinement_->refine(1, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                     valueOf(getProblem()->getSolution(0))));*/
    /*								 
	if(!(rank_ == 2 && time_counter_<10)){
		refinement_->refine(1, function_(indicator(getMesh()->getName() + "->phase", 3 * eps_),
                                     valueOf(getProblem()->getSolution(0))));
	}*/
    MSG("time (refinement) = %e\n", t.elapsed());

    //rescaling phase field
    //*getProblem()->getSolution(0) << max(-1.0, min(1.0, valueOf(*getProblem()->getSolution(0))));


    // Possibly write the velocity ...
    bool writeVelocities = true;
    Parameters::get("main->write velocities", writeVelocities);

    // Update cell position
    oldPos_ = pos_;
    pos_ = integrate(X() * ((valueOf(*tmp_lagrange2_) + 1.0) / 2.0)) / integrate((valueOf(*tmp_lagrange2_) + 1.0) / 2.0);
    
    updateVelocity();

    if (writeVelocities)
    {
      WorldVector<double> velocity = (pos_- oldPos_) * (1.0 / adaptInfo->getTimestep());

      std::string filename = directory_ + "/velocities" + postfix_ + ".dat";
      std::ofstream velocityFile;
      velocityFile.open(filename, std::ios::app);
      velocityFile << adaptInfo->getTime() << " " << velocity << "\n";

      velocityFile.close();
    }

    if (writeData_){   //maybe this line is changed
      std::ofstream out(filename_positions_, std::ios_base::app);
      //WorldVector<double> S;
      //S[0] = integrate(valueOf(nematic_tensor_[0]));
      //S[1] = integrate(valueOf(nematic_tensor_[1]));
      bool center_of_mass = 1;
      Parameters::get("main->center of mass", center_of_mass);
      WorldVector<double> center;

      if (center_of_mass)
          center = getPosition();
      else
          center = getCenter();
      
      if(confined_ == 1){
        confinement_interaction_ = integrate(valueOf(*confinementForce_))/Con_;
      }
      out << adaptInfo->getTime() << ',' << rank_ << ',' << center[0] << ',' << center[1] 
        << ',' << radius_ << ',' << nematic_tensor_[0] << ',' << nematic_tensor_[1] << ',' 
        << velocityValues_[0] << ',' << velocityValues_[1] << ',' << major_axis_angle_ 
        << ','<< total_interactions << ',' << neighbours << ','<< confinement_interaction_ << ',' << actual_growth_factor_ << '\n';
    }
    Super::closeTimestep(adaptInfo);
  }

  // After the position was updated we can compute the velocity
  void updateVelocity()
  {
    WorldVector<double> vel;
    vel[0] = pos_[0] - oldPos_[0];
    vel[1] = pos_[1] - oldPos_[1];

    if (std::abs(vel[0]) > 0.5 * domainDimension_[0])
    { 
      vel[0] = vel[0]>0?vel[0]-domainDimension_[0]:vel[0]+domainDimension_[0];
    }
    if (std::abs(vel[1]) > 0.5 * domainDimension_[1])
    {
      vel[1] = vel[1]>0?vel[1]-domainDimension_[1]:vel[1]+domainDimension_[1];
    }
    velocityValues_[0] = vel[0];
    velocityValues_[1] = vel[1];
  }

  // expressiong representing the Double-Well
  template <class Phi>
  auto B(Phi &&p)
  {
    return max(0.01,  9.0/4.0*pow<2>(pow<2>(std::forward<Phi>(p)) - 1.0));
    //return 1.0;
  }
  template <class Phi>

  auto B2(Phi &&p)
  { 
  return max(0.00, pow<2>(pow<2>(std::forward<Phi>(p)) - 1.0));
  }

  //Stabilization
  template <class Phi>
  auto G(Phi &&p)
  { 
  return sqrt(4.0/9.0*pow<2>(std::forward<Phi>(p)-1.0)*pow<2>(std::forward<Phi>(p)+1.0));
  //return sqrt(4.0/9.0*pow<2>(std::forward<Phi>(p)-1.0)*pow<2>(std::forward<Phi>(p)+1.0)+eps_*eps_*1.0e-1);
  }


  virtual void fillOperators(ProblemType *prob) override
  {
	
    auto *phi = prob->getSolution(0);
    auto phi_restrict = max(-1.0, min(1.0, valueOf(phi)));
    //auto *growth_f = prob.getGrowthFactor();
    auto *advec = &advection_;
    // dt(phi) - gamma*laplace(phi^#) = -lambda1
    // -----------------------------------------

    prob->addTimeOperator(0, 0, getInvTau());
	/*
    Operator *T = new Operator(prob->getFeSpace(0), prob->getFeSpace(0));
    T->setUhOld(*getProblem->getOldSolution(0));
    addZOT(T,1.0);
    prob->addMatrixOperator(T,0,0, getInvTau(), getInvTau());
    prob->addVectorOperator(T,0, getInvTau());
	*/

    Operator *opLaplace = new Operator(prob->getFeSpace(0), prob->getFeSpace(1));
    addSOT(opLaplace,gamma_);///conservative B(valueOf(phi))*B(valueOf(phi))*
    //addZOT(opLaplace, gamma2_);///non conservative
    //addZOT(opLaplace, allen_growth_);
    prob->addMatrixOperator(opLaplace, 0, 1);

    // Random advection term
    Operator *opAdvection = new Operator(prob->getFeSpace(0));
    //addFOT(opAdvection, valueOf(advection_), GRD_PSI);
    //addFOT(opAdvection, valueOf(adRef()), GRD_PSI);//e
    //addFOT(opAdvection, valueOf(adref_), GRD_PSI);//f
    addFOT(opAdvection, valueOf(*advec), GRD_PSI);

    prob->addMatrixOperator(opAdvection, 0, 0);  

    // phi^# = 1/eps * (phi^3 - phi) - eps*laplace(phi)
    // ------------------------------------------------

    Operator *opPhi2 = new Operator(prob->getFeSpace(1));
    addZOT(opPhi2, G(valueOf(phi))*1.0);//G(valueOf(phi))*G(phi_restrict)*G(valueOf(phi))*

    prob->addMatrixOperator(opPhi2, 1, 1);

    Operator *opF_lhs = new Operator(prob->getFeSpace(1), prob->getFeSpace(0));
    addZOT(opF_lhs, (1.0 / (Ca_* eps_)) * (1.0 - 3.0 * pow<2>(valueOf(phi))));

    Operator *opF_rhs = new Operator(prob->getFeSpace(1));
    addZOT(opF_rhs, (-2.0 / (Ca_ * eps_)) * pow<3>(valueOf(phi)));

    prob->addMatrixOperator(opF_lhs, 1, 0);
    prob->addVectorOperator(opF_rhs, 1);

    Operator *opLaplace2 = new Operator(prob->getFeSpace(1), prob->getFeSpace(0));
    addSOT(opLaplace2, -eps_ / (Ca_));

    prob->addMatrixOperator(opLaplace2, 1, 0);

    //operator for growth
    Operator *opGrowth = new Operator(prob->getFeSpace(0));
    addZOT(opGrowth, ref_(f_growth_)*(valueOf(phi)+1)/2.0);//ref_(f_growth_)(*growth_f)
    prob->addVectorOperator(opGrowth, 0);
    
    //Operator *opGrowthM = new Operator(prob->getFeSpace(0), prob->getFeSpace(0));
    //addZOT(opGrowthM, -ref_(f_growth_));
    //prob->addMatrixOperator(opGrowthM, 0, 0);
  

    //Operator *opSource = new Operator(prob->getFeSpace(0));
    //addZOT(opSource, So_);
    //prob->addVectorOperator(opSource, 0);

    if (confined_==1){
      // B'(phi) * w(d_wall)
      Operator *opConf_rhs = new Operator(prob->getFeSpace(1));
      addZOT(opConf_rhs, 1.0/Con_*valueOf(*confinementForce_)*G(valueOf(phi)));
      prob->addVectorOperator(opConf_rhs, 1);
    }

  }

  using Super::fillOperators;

  virtual void fillBoundaryConditions(ProblemType *prob) override
  {
    // priodic boundary conditions for left-right (-1) and bottom-top (-2) periodicity
    for (int j = 0; j < prob->getNumComponents(); ++j)
    {
      for (int k = 0; k < prob->getNumComponents(); ++k)
      {
        prob->addPeriodicBC(-1, j, k);
        prob->addPeriodicBC(-2, j, k);
      }
    }
    //prob->addDirichletBC(1, 0, 0, new G);
    
  }

  using Super::fillBoundaryConditions;

  auto const &getVerletList() const { return verletList_->get(); }
  
  //auto const &getVerletListNoOverlap() const { return verletList_->getNoOverlap(); }

  // Computed position, i.e. center of mass
  WorldVector<double> getPosition() const
  {
    return pos_;  
  }
  
  // Position of central DOF
  WorldVector<double> getCenter() const
  {
    return verletList_->getCenter();
  }

  // return upper bound of radius of cell
  double getRadius() const
  {
    return radius_;
  }

  double getVolume() const
  {
    return volume_;
  }

  WorldVector<DOFVector<double>*> adRef() const{
    return advection_;
  }

  auto adRef2() {
    return &(advection_);
  }
  /*
  WorldVector<DOFVector<double>*> adChange(){
    return valueOf(advection_[0])*ref_(ct), valueOf(advection_[1])*ref_(st);
  }*/

  WorldVector<double> getNematicTensor() const
  {
    return nematic_tensor_;
  }

  double* getGrowthFactor(){
    return &f_growth_;
  }



  class G : public AbstractFunction<double, WorldVector<double> >
  {
  public:
    double operator()(const WorldVector<double>& x) const
    {
      return -1.0;
    }
  };

protected:
  std::unique_ptr<DOFVector<double>> tmp_lagrange2_, tmp_lagrange1_;

  std::unique_ptr<VerletList> verletList_;
  std::unique_ptr<MacroCommunicator> communicator_;
  
  double radius_=1.0;

  double initial_volume_, volume_, current_volume_;

  WorldVector<double> pos_,oldPos_;
  WorldVector<double> domainDimension_;
  std::unique_ptr<HL_SignedDistTraverse> reinit_;
  std::unique_ptr<extensions::RefinementExpression> refinement_;

  std::unique_ptr<ContourWriter> writer_signedDist_;

  double h_ = 1.0; // grid-width

  double gamma_ = 1.0;
  double eps_ = 0.2;
  double In_ = 0.1125;
  double Ca_ = 0.0281;
  double rescaledCa_ = Ca_;
  double rescaledIn_ = In_;
  int div_refine_ = 0;
  
  int neighbours = 0; //Number of cells, a cell is interacting with

  double v0_ = 0.5;
  double allen_growth_ = 0.0;
  double mu_laplace_const_ = 1.0;
  double volume_threshold_ = 1.1e3;
  bool boost_refinement_ = false;
  double growth_factor_ = 0.1;
  double actual_growth_factor_ = 0.1;
  double confinement_inhibition_ = 1.0; //1 = no inhibition, 0.0 = complete inhibition
  double con_inhibiting_limit_ = 500.0;
  int growth_type_ = 1;
  double f_growth_ = 1.0;
  double &f_growth_ref{f_growth_};
  double Gr_ = 1.0;
  double gamma2_ = 1.0;
  double So_;
  int refinement_counter = 0;

  double total_interactions = 0.0;

  double resize = 0.25;

  int time_threshold_ = 10;
  int time_counter_ = 0;
  int numInitialCells_ = 1;
  int phaseStatus = 0;
  std::vector<int> rankTree_;
  int parent_;
  int nextChildIndex_;

  std::vector<int> aliveCells_; //list of all cells which are alive
  std::vector<int> waitingCells_; //list of cells waiting to be born
  int latestWaitingCell_;
  int alive_ = 0; //1 = alive, 2 = dead, 0 = not born yet

  bool writeData_ = false;

  std::string directory_, postfix_, setup_, filename_positions_;

  // nematic deformation tensor
  WorldVector<double> nematic_tensor_;
  double major_axis_angle_;

  // Advection stuff
  //std::unique_ptr<WorldVector<DOFVector<double>>> advection_;
  WorldVector<DOFVector<double>*> advection_;
  WorldVector<DOFVector<double>*>& adref_ = advection_; 
  double theta_;
  double angleDiffusivity_; // How much can we change in a single step?
  
  // velocity
  WorldVector<double> velocityValues_;

  //confinement
  std::unique_ptr<DOFVector<double>> confinementPF_; //confinement phasefield
  std::unique_ptr<DOFVector<double>> confinementForce_;//a DOF Vector that stores repulsive contributions from confinement for each cell
  int confinementShape_ = 1; //1: rectangle or 2: circle 
  int confined_ = 0;
  double Con_;
  double confinement_interaction_ = 0.0;
  double confinement_size_ = 1.0;

  double timesteptest = 0;

  int random_seed = 1020;
  double ct = 0.0;
  double st = 0.0;

  //variable to set size of one of the two cells to zero
  int switchsize = 1; //if 1 then cell 1 exists, 2 then cell 2 exists and 3 then both exist
};

} // namespace base_problems
} // namespace AMDiS
