#ifndef AMDIS_ADAPT_INSTATIONARY_BASE_H
#define AMDIS_ADAPT_INSTATIONARY_BASE_H

#include "AMDiS.h"

#include <mpi14/mpi14.hpp>

namespace AMDiS {

  class AdaptInstationaryBase : public AdaptInstationary
  {
  public:
    AdaptInstationaryBase(const std::string name,
			ProblemIterationInterface &problemStat,  
			AdaptInfo &adaptInfo,
			ProblemTimeInterface &problemInstat,
			AdaptInfo &initialInfo,
			std::time_t initialTimestamp=0) 
      :	AdaptInstationary(name, problemStat, adaptInfo, problemInstat, initialInfo, initialTimestamp),
	breakConditionTemp(false),
	minTimesteps(0)
    {
      breakCondition= &breakConditionTemp;
      Initfile::get("user parameter->min timesteps", minTimesteps, 2);
    }
	    
    virtual void reset(AdaptInfo *adaptInfo) 
    {
      adaptInfo->reset();
      adaptInfo->init();
      *breakCondition= false;
    }
    
    virtual void writeSolution(int problemNr) {}
    virtual void initAdaption(AdaptInfo* adaptInfo) {}
    
    void setBreakCondition(bool *breakCondition_) 
    {
      breakCondition=breakCondition_;
    }
    
    int adapt()
    {
	    
      int errorCode = 0;
      
      TEST_EXIT(adaptInfo->getTimestep() >= adaptInfo->getMinTimestep())
	("timestep < min timestep\n");
      TEST_EXIT(adaptInfo->getTimestep() <= adaptInfo->getMaxTimestep())
	("timestep > max timestep\n");
	    
      #if HAVE_PARALLEL_DOMAIN_AMDIS
      Parallel::MeshDistributor::globalMeshDistributor->initParallelization(); 
      #endif

      if (adaptInfo->getTimestepNumber() == 0) {
	adaptInfo->setTime(adaptInfo->getStartTime());
	initialAdaptInfo->setStartTime(adaptInfo->getStartTime());
	initialAdaptInfo->setTime(adaptInfo->getStartTime());

	problemTime->setTime(adaptInfo);

	// initial adaption
	problemTime->solveInitialProblem(initialAdaptInfo);
	problemTime->transferInitialSolution(adaptInfo);

	initAdaption(adaptInfo);
      }
      
      if (adaptInfo->getStartTime() >= adaptInfo->getEndTime()) {
	problemTime->closeTimestep(adaptInfo);
	return 0;
      }
	    
      while (!adaptInfo->reachedEndTime()) {
	iterationTimestamp = time(NULL);

        MSG("start (AdaptInstationary::initTimestep) = %20.16e\n", mpi14::now());
	problemTime->initTimestep(adaptInfo);
        
        MSG("start (AdaptInstationary::oneTimestep) = %20.16e\n", mpi14::now());
	oneTimestep();
        
        MSG("start (AdaptInstationary::closeTimestep) = %20.16e\n", mpi14::now());
	problemTime->closeTimestep(adaptInfo);

	if (breakWhenStable && (adaptInfo->getSolverIterations() == 0))
	  break;
	
	if(*breakCondition && adaptInfo->getTimestepNumber() > minTimesteps) {
	  adaptInfo->setTime(adaptInfo->getEndTime());
	  adaptInfo->setTimestepNumber(adaptInfo->getNumberOfTimesteps());
	  break;
	}

	// Check if there is a runtime limitation. If there is a runtime limitation
	// and there is no more time for a next adaption loop, than return the error
	// code for rescheduling the problem and break the adaption loop.
	if (checkQueueRuntime()) {
	  errorCode = RescheduleErrorCode;
	  break;
	}
      }
	    
      return errorCode;
    }
	  
  protected:
    bool *breakCondition;
    bool breakConditionTemp;
    int minTimesteps;
  };

} // end namespace AMDiS

#endif // AMDIS_ADAPT_INSTATIONARY_BASE_H
