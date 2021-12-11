#pragma once

#include <string>

#include <mpi14/mpi14.hpp>

#include "AdaptInfo.h"
#include "ProblemIterationInterface.h"
#include "ProblemInstat.h"

namespace AMDiS {
  namespace base_problems {

/**
 * An MPIBaseProblem is a ScalarBaseProblem, that has an additional identifier, namely an process-id
 **/
template <class ProblemType>
class MPIBaseProblem
    : public ProblemIterationInterface
    , public ProblemInstatBase
{
public:

  MPIBaseProblem(std::string const& name, mpi14::Communicator comm = {})
    : ProblemInstatBase(name, nullptr)
    , name_(name)
    , problem_(name + "->space")
    , comm_(comm)
  {
    rank_ = comm_.rank();
    size_ = comm_.size();

    dow_ = Global::getGeo(WORLD);
    Parameters::get(name_ + "->space->dim", dim_);

    appendOutputPostfix("p" + std::to_string(rank_));
  }

  /// Initialisation of the problem.
  virtual void initialize(Flag initFlag,
                          ProblemStat *adoptProblem = NULL,
                          Flag adoptFlag = INIT_NOTHING)
  {
    problem_.initialize(initFlag, adoptProblem, adoptFlag);
  }

  /// Initialisation of DOFVectors and AbstractFunctions,
  /// is called in \ref initTimeInteface after feSpace and mesh are initialized
  virtual void initData()
  {}

  /// Method is called at the end of \ref initTimeInteface
  virtual void finalizeData() {}

  /// calls \ref initData, \ref fillOperators and \ref fillBoundaryConditions in this ordering
  virtual void initBaseProblem()
  {
    initData();
    fillOperators();
    fillBoundaryConditions();
    finalizeData();
  }

  virtual void solveInitialProblem(AdaptInfo *adaptInfo) override
  {}

  /// calls \ref writeFiles and updates \ref oldMeshChangeIdx_
  virtual void transferInitialSolution(AdaptInfo *adaptInfo) override
  {
    writeFiles(adaptInfo, false);
  }

  /// This method is called before \ref beginIteration, \ref oneIteration and \ref endIteration.
  virtual void initTimestep(AdaptInfo *adaptInfo) override
  {}

  /// calls \ref writeFiles
  virtual void closeTimestep(AdaptInfo *adaptInfo) override
  {
    writeFiles(adaptInfo, false);
  }

  virtual void beginIteration(AdaptInfo *adaptInfo) override
  { FUNCNAME("MPIBaseProblem::beginIteration()");

    MSG("\n");
    MSG(("[[ <"+name_+"> iteration ]]\n").c_str());
  }

  virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo = FULL_ITERATION) override
  {
    Flag flag;

    problem_.beginIteration(adaptInfo);
    flag |= problem_.oneIteration(adaptInfo, toDo);
    problem_.endIteration(adaptInfo);

    return flag;
  }

  virtual void endIteration(AdaptInfo *adaptInfo) override
  { FUNCNAME("MPIBaseProblem::endIteration()");

    MSG("\n");
    MSG(("[[ end of <"+name_+"> iteration ]]\n").c_str());
  }

  /// Calls writeFiles of the problem
  virtual void writeFiles(AdaptInfo *adaptInfo, bool force)
  {
    problem_.writeFiles(adaptInfo, force);
  }

  // getting methods

  /// pointer to the mesh of the problem
  Mesh* getMesh(int comp = 0)
  {
    return problem_.getMesh(comp);
  }

  Mesh const* getMesh(int comp = 0) const
  {
    return problem_.getMesh(comp);
  }

  /// pointer to the feSpace of the problem
  const FiniteElemSpace* getFeSpace(int comp = 0) const
  {
    return problem_.getFeSpace(comp);
  }

  /// name of the baseBroblem
  std::string getName() const
  {
    return name_;
  }

  int getRank() const 
  {
    return rank_;
  }

  int getNumProblems() const
  {
    return 1;
  }

  int getNumComponents() const
  {
    return problem_.getNumComponents();
  }

  ProblemType* getProblem(int=0)
  {
    return &problem_;
  }

  ProblemType const* getProblem(int=0) const
  {
    return &problem_;
  }

  // setting methods
  void serialize(std::ostream&) {}

  void deserialize(std::istream&) {}

  virtual void fillOperators(ProblemType* /*prob*/)
  {}

  /// method where operators are added to the problem
  virtual void fillOperators()
  {
    fillOperators(&problem_);
  }

  virtual void fillBoundaryConditions(ProblemType* /*prob*/)
  {}

  /// method where boundary conditions are added to the problem
  virtual void fillBoundaryConditions()
  {
    fillBoundaryConditions(&problem_);
  }

  mpi14::Communicator const& comm() const { return comm_; }
  mpi14::Communicator& comm() { return comm_; }

protected:

  void appendOutputPostfix(std::string postfix) const
  {
    // copy parameter values from main problem to individual problems
    std::map<std::string,std::string> params;

    std::string base_name = name_ + "->space";
    Parameters::getParameterMap(base_name + "->", params);

    if (params.count("output->filename") > 0) {
      std::string filename_i = params["output->filename"] + postfix + "_";
      Parameters::set(base_name + "->output->filename", filename_i);
    }

    int components = 0;
    Parameters::get(base_name + "->components", components);
    for (int c = 0; c < components; ++c) {
      std::string output_i = "output[" + std::to_string(c) + "]";
      if (params.count(output_i + "->filename") > 0) {
        std::string filename_i = params[output_i + "->filename"] + postfix + "_";
        Parameters::set(base_name + "->" + output_i + "->filename", filename_i);
      }
    }

    int vectors = 0;
    Parameters::get(base_name + "->output->num vectors", vectors);
    for (int v = 0; v < vectors; ++v) {
      std::string output_i = "output->vector[" + std::to_string(v) + "]";
      if (params.count(output_i + "->filename") > 0) {
        std::string filename_i = params[output_i + "->filename"] + postfix + "_";
        Parameters::set(base_name + "->" + output_i + "->filename", filename_i);
      }
    }
  }

protected:

  /// The name of this base-problem, used in init-files
  std::string name_;

  /// The underlying problemStat
  ProblemType problem_;

  /// The MPI communicator to use, e.g. MPI_COMM_WORLD
  mpi14::Communicator comm_;

  int rank_ = 0; // processor id
  int size_ = 1; // num processors

  /// dimension of the mesh (set in \ref initialize(...) )
  int dim_ = -1;

  /// dimension of world
  int dow_ = Global::getGeo(WORLD);

public:
  /// additions
  int nextEmptyRank = 1;
  int latestParentRank = 0;
};
} }
