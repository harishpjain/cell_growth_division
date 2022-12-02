#pragma once

#include "AMDiS_fwd.h"
#include "AdaptInstationaryBase.h"

#include <ctime>
#include <string>
#include <vector>

namespace AMDiS 
{

  class AdaptInstationarySeq 
      : public AdaptInstationaryBase
  {
  public:
    AdaptInstationarySeq(const std::string& name,
			ProblemIterationInterface& problemStat,  
			AdaptInfo& adaptInfo,
			ProblemTimeInterface& problemInstat,
			AdaptInfo& initialAdaptInfo,
			std::time_t initialTimestamp = 0);
			
    void initAdaption(AdaptInfo* adaptInfo);
    
    void setNewTimestep();
    
    void timeStatistics(Timer& t);
    
    void explicitTimeStrategy();
    
    void implicitTimeStrategy();
    
    int getNrOfIterations() { return totalNrOfTimesteps_; }
    
    void reset(AdaptInfo *adaptInfo);

    void writeSolution(int problemNr);

  private:
    bool useEstimator_ = false;
    
    int totalNrOfTimesteps_ = 0;
    int adaptAfterTimesteps_ = 0;
    int info_ = 10;
    
    double totalTime_ = 0.0;
    double avIterTime_ = 0.0;
    double endTime_ = 0.0;
    
    std::vector<double> timesteps_;
    std::vector<int> nrOfTimesteps_;
  };

} // end namespace AMDiS
