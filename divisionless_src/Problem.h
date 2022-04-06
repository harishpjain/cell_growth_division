#pragma once

#include "AMDiS.h"
#include "ProblemStat.h"
#include "extensions/RefineOperations.h"

namespace AMDiS {

  class Problem
      : public ProblemStat
      , public extensions::StandardRefineOperation
  {
  public:
    Problem(std::string name_)
      : ProblemStat(name_)
    {}

    void addTimeOperator(int i, int j, double* invTau)
    {
      oldSolution_.reset(new DOFVector<double>(this->getFeSpace(j), "old_solution_" + std::to_string(j)));

      Operator* opTime = new Operator(this->getFeSpace(i), this->getFeSpace(j));
      addZOT(opTime, 1.0);
      opTime->setUhOld(oldSolution_.get());
      this->addMatrixOperator(opTime, i,j, invTau);
      this->addVectorOperator(opTime, i, invTau);
    }

    DOFVector<double>* getOldSolution() { return oldSolution_.get(); }

  private:
    std::unique_ptr<DOFVector<double>> oldSolution_;
  };

} // end namespace AMDiS
