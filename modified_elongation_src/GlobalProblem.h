#pragma once

#include "MPIBaseProblem.h"

#include "extensions/RefinementExpression.h"
#include "extensions/Views.h"

namespace AMDiS { namespace base_problems {

class GlobalProblem
    : public MPIBaseProblem<ProblemStat>
{
  using Super = MPIBaseProblem<ProblemStat>;

public:
  GlobalProblem(std::string name, mpi14::Communicator comm)
    : Super(name, comm)
  {
    int comp = 0;
    Parameters::get(name + "->space->components", comp);
    for (int c = 0; c < comp; ++c) {
      int ref = -1 - c;
      Parameters::set(name + "->space->refinement set[" + std::to_string(c) + "]", ref);
    }
  }

  virtual void beginIteration(AdaptInfo */*adaptInfo*/) override {}
  virtual Flag oneIteration(AdaptInfo */*adaptInfo*/, Flag /*toDo*/ = FULL_ITERATION) override
  {
    Flag null;
    return null;
  }
  virtual void endIteration(AdaptInfo */*adaptInfo*/) override {}
};

}} // end namespaces
