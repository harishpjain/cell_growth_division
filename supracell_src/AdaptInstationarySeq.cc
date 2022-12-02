#include <algorithm>
#include <numeric>

#include "AdaptInstationarySeq.h"
#include "Helpers.h"

namespace AMDiS {

AdaptInstationarySeq::AdaptInstationarySeq(const std::string& name,
				ProblemIterationInterface& problemStat,  
				AdaptInfo& info,
				ProblemTimeInterface& problemInstat,
				AdaptInfo& initialInfo,
				std::time_t initialTimestamp) 
  : AdaptInstationaryBase(name,problemStat,info,problemInstat,initialInfo,initialTimestamp)
{
  FUNCNAME("AdaptInstationarySeq::AdaptInstationarySeq()");
  Parameters::get("user parameter->use estimator", useEstimator_, 2);

  Parameters::get(name + "->sequence->timesteps", timesteps_, 2);
  Parameters::get(name + "->sequence->number of timesteps", nrOfTimesteps_, 2);
  Parameters::get(name + "->start after first timesteps", adaptAfterTimesteps_, 2);

  TEST_EXIT(timesteps_.size() == nrOfTimesteps_.size() && timesteps_.size() > 0)
    ("#('timesteps') and #('number of timesteps') must be equal! (%d /= %d)\n", 
     timesteps_.size(), nrOfTimesteps_.size());

  std::string directory = ".";
  std::string postfix = "";
  Parameters::get("output", directory);
  Parameters::get("postfix", postfix);
  
  std::remove((directory + "/timeStatistics" + postfix + ".csv").c_str());
}

void AdaptInstationarySeq::initAdaption(AdaptInfo* adaptInfo)
{ FUNCNAME("AdaptInstationarySeq::initAdaption()");

  // set endTime for simulation
  endTime_ = std::inner_product(timesteps_.begin(), timesteps_.end(), nrOfTimesteps_.begin(), 0.0);

  bool calcEndTime = false;
  Parameters::get(name+"->sequence->calc end time", calcEndTime, 2);
  if (!calcEndTime && adaptInfo->getEndTime() < endTime_)
    endTime_ = adaptInfo->getEndTime();
  adaptInfo->setEndTime(endTime_);

  adaptInfo->setMaxTimestep(*std::max_element(timesteps_.begin(), timesteps_.end()));
  totalNrOfTimesteps_ = std::accumulate(nrOfTimesteps_.begin(), nrOfTimesteps_.end(), 0);
  if (adaptInfo->getStartTime() >= adaptInfo->getEndTime())
    totalNrOfTimesteps_ = 0;
  
  adaptInfo->setNumberOfTimesteps(totalNrOfTimesteps_);

  INFO(2,2)("===========================\n");
  INFO(2,2)("EndTime: %f\n", endTime_);
  INFO(2,2)("NrOfIterations <= %d\n", totalNrOfTimesteps_);
  INFO(2,2)("===========================\n");
}

void AdaptInstationarySeq::setNewTimestep()
{ FUNCNAME("AdaptInstationarySeq::setNewTimestep()");

  if (totalNrOfTimesteps_ <= 0)
    return;
  
  std::vector<int> barriers(nrOfTimesteps_.size());
  std::partial_sum(nrOfTimesteps_.begin(), nrOfTimesteps_.end(), barriers.begin());

  double newTimestep = 1.0;
  for (std::size_t i = 0; i < barriers.size(); ++i) {
    if (adaptInfo->getTimestepNumber() < barriers[i]) {
      newTimestep = timesteps_[i];
      break;
    }
  }
  
  if (adaptInfo->getTime() + newTimestep >= adaptInfo->getEndTime() - DBL_TOL) {
    adaptInfo->setNumberOfTimesteps(adaptInfo->getTimestepNumber());
    adaptInfo->setEndTime(adaptInfo->getTime() + newTimestep);
  }
  adaptInfo->setTimestep(newTimestep); 

  if (newTimestep < std::max(1.e-10, adaptInfo->getMinTimestep())) {
    WARNING("last timestep <= %e. Stop before endTime!\n", adaptInfo->getMinTimestep());
    adaptInfo->setEndTime(adaptInfo->getTime());
    return;
  }
  
  adaptInfo->setTime(adaptInfo->getTime() + newTimestep);	
  INFO(2,1)("time = %e, timestep = %e\n", adaptInfo->getTime(), newTimestep);
  problemTime->setTime(adaptInfo);
}

void AdaptInstationarySeq::timeStatistics(Timer& t)
{ FUNCNAME("AdaptInstationarySeq::timeStatistics()");

  // statistic about time left
  double iterTime = t.elapsed();
  int iter = adaptInfo->getTimestepNumber();
  totalTime_ += iterTime;
  avIterTime_ = totalTime_/(iter+1);

  MSG("\n");
  MSG(("~~~~~(%d/%d)" + Helpers::fillString(42,'~',2,iter+1,totalNrOfTimesteps_) + "\n").c_str(),iter+1,totalNrOfTimesteps_);
  MSG("iter-time: %.5fs, [average: %.5fs]\n",iterTime,avIterTime_);
  MSG("remaining time: %.5f min\n",((totalNrOfTimesteps_-iter-1)*avIterTime_/60.0));
  MSG("#iterations left: %d\n",(totalNrOfTimesteps_-iter-1));
  MSG((Helpers::fillString(50,'~',0) + "\n").c_str());
}

void AdaptInstationarySeq::explicitTimeStrategy()
{
  setNewTimestep();
  adaptInfo->setSpaceIteration(0);

  // do the iteration
  Timer t;
  problemIteration->beginIteration(adaptInfo);
  Flag iterFlag = FULL_ITERATION;
  if (!useEstimator_) {
    iterFlag.unsetFlag(ADAPT);
    iterFlag.unsetFlag(MARK);
    iterFlag.unsetFlag(ESTIMATE);
  }
  problemIteration->oneIteration(adaptInfo, iterFlag);
  problemIteration->endIteration(adaptInfo);

  timeStatistics(t);
  adaptInfo->setLastProcessedTimestep(adaptInfo->getTimestep());
}

void AdaptInstationarySeq::implicitTimeStrategy()
{
  setNewTimestep();
  adaptInfo->setSpaceIteration(0);

  Timer t;

  problemIteration->beginIteration(adaptInfo);
  problemIteration->oneIteration(adaptInfo, NO_ADAPTION);
  // Do only space iterations only if the maximum is higher than 0
  if (adaptInfo->getMaxSpaceIteration() > 0 &&
      adaptInfo->getTimestepNumber() >= adaptAfterTimesteps_) {
    do { // Space iterations
      problemIteration->oneIteration(adaptInfo, FULL_ITERATION);
      adaptInfo->incSpaceIteration();
    } while (!adaptInfo->spaceToleranceReached() &&
	      adaptInfo->getSpaceIteration() <= adaptInfo->getMaxSpaceIteration());
  }
  problemIteration->endIteration(adaptInfo);

  timeStatistics(t);
  adaptInfo->setLastProcessedTimestep(adaptInfo->getTimestep());
}

void AdaptInstationarySeq::reset(AdaptInfo *adaptInfo)
{
  AdaptInstationaryBase::reset(adaptInfo);
  adaptInfo->setTimestep(timesteps_[0]);
  adaptInfo->setTimestepNumber(0);
  totalTime_ = 0.0;
  avIterTime_ = 0.0;
}

void AdaptInstationarySeq::writeSolution(int problemNr) 
{
  if (problemNr > 0 && problemNr < problemIteration->getNumProblems())
    static_cast<ProblemStat*>(problemIteration->getProblem(problemNr))->writeFiles(adaptInfo, true);
}

} // end namespace AMDiS
