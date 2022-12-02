/** \file CouplingBaseProblem.h */

#ifndef COUPLING_BASE_PROBLEM_H
#define COUPLING_BASE_PROBLEM_H

#include "AMDiS.h"

// coupling structures
#include "CouplingIterationInterface.h"
#include "CouplingTimeInterface.h"
#include "CouplingProblemStat.h"
#include "BaseProblem.h"
#include "ProblemInstat.h"

#include <tuple>
#include "GenericLoops_cxx11.h"

namespace AMDiS {

namespace extensions {

  namespace detail {

    /// for generic loops this struct can be passed
    struct AddProblem {
      template<typename CouplingProblemStatType, typename BaseProblemType>
      static void call(CouplingProblemStatType& couplingProblemStat, BaseProblemType& baseProblem) {
	for (size_t j = 0; j < baseProblem.getNumProblems(); j++)
	  couplingProblemStat.addProblem(baseProblem.getProblem(j));
      }
    };

    /// for generic loops this struct can be passed
    struct AddTimeInterface {
      template<typename CouplingProblemStatType, typename BaseProblemType>
      static void call(CouplingProblemStatType& couplingProblemStat, BaseProblemType& baseProblem) {
	couplingProblemStat.addTimeInterface(&baseProblem);
      }
    };

    /// for generic loops this struct can be passed
    struct AddIterationInterface {
      template<typename CouplingProblemStatType, typename BaseProblemType>
      static void call(CouplingProblemStatType& couplingProblemStat, BaseProblemType& baseProblem) {
	couplingProblemStat.addIterationInterface(&baseProblem);
      }
    };

    /// Functor for generic loops. Method initData() is called for each element in a sequence.
    struct InitData {
      template<typename BaseProblemType>
      static void call(BaseProblemType& baseProblem) { baseProblem.initData(); }
    };

    /// Functor for generic loops. Method initData() is called for each element in a sequence.
    struct FinalizeData {
      template<typename BaseProblemType>
      static void call(BaseProblemType& baseProblem) { baseProblem.finalizeData(); }
    };

    /// Functor for generic loops. Method fillOperators() is called for each element in a sequence.
    struct FillOperators {
      template<typename BaseProblemType>
      static void call(BaseProblemType& baseProblem) { baseProblem.fillOperators(); }
    };

    /// Functor for generic loops. Method fillBoundaryConditions() is called for each element in a sequence.
    struct FillBoundaryConditions {
      template<typename BaseProblemType>
      static void call(BaseProblemType& baseProblem) { baseProblem.fillBoundaryConditions(); }
    };


    struct FindProblem {
      template<typename BaseProblemType, typename ProblemType>
      static void call(BaseProblemType& baseProblem, const std::string& name, ProblemType& prob) {
	if (baseProblem.getName() == name)
	  prob = baseProblem.getProblem();
      }
    };

    struct FindBaseProblem {
      template<typename BaseProblemType, typename ProblemType>
      static void call(BaseProblemType& baseProblem, const std::string& name, ProblemType*& prob) {
	typedef typename boost::mpl::if_<typename boost::is_same<BaseProblemType, ProblemType>::type,
					boost::mpl::bool_<true>,
					boost::mpl::bool_<false> >::type assign;
	call(baseProblem, name, prob,  assign());
      }
      template<typename BaseProblemType, typename ProblemType>
      static void call(BaseProblemType& baseProblem, const std::string& name, ProblemType*& prob, boost::mpl::bool_<true>) {
	if (baseProblem.getName() == name)
	  prob = &baseProblem;
      }
      template<typename BaseProblemType, typename ProblemType>
      static void call(BaseProblemType& baseProblem, const std::string& name, ProblemType*& prob, boost::mpl::bool_<false>) {}
    };

  } // end namespace detail

/**
  * \ingroup Problem
  *
  * \brief Structur to couple BaseProblems of variouse types
  */
template <typename ProblemType=ProblemStat, typename... BaseProblemTypes>
class CouplingBaseProblem
    : public CouplingIterationInterface
    , public CouplingTimeInterface
    , public AMDiS::detail::CouplingProblemStat<ProblemType>
{
  typedef AMDiS::detail::CouplingProblemStat<ProblemType> CProblemStat;
  typedef std::tuple<BaseProblemTypes&...> BaseProblemsTupleType;

public:
  CouplingBaseProblem(std::string name_, BaseProblemTypes&... baseProblems_)
    : ProblemIterationInterface()	// virtual base class constructor
    , ProblemTimeInterface()		// virtual base class constructor
    , CProblemStat(name_)
    , baseProblems(baseProblems_...)
    , name(name_)
  {
    dow = Global::getGeo(WORLD);
    Parameters::get(name_ + "->dim", dim);
    
    tools::FOR_EACH< detail::AddProblem >::loop2(*this, baseProblems);
    tools::FOR_EACH< detail::AddIterationInterface >::loop2(*this, baseProblems);
    tools::FOR_EACH< detail::AddTimeInterface >::loop2(*this, baseProblems);
  }

  virtual ~CouplingBaseProblem() { }

  /**
   * Add the problems to the iterationInterface, timeInterface and couplingProblemStat.
   * As a consequence all problem can be initialized one after another and in the
   * adaption loop they are solved in rotation.
   *
   * In the adaption loop the problems are solved the same order as they are added to the
   * iterationInterface in this method. This order can be changed manually in the oneIteration
   * method.
   **/
  virtual void initialize(Flag initFlag,
		  ProblemStatSeq *adoptProblem = NULL,
		  Flag adoptFlag = INIT_NOTHING) override
  {    
    // initialize all ProblemStat
    CProblemStat::initialize(initFlag, adoptProblem, adoptFlag);
  }

  virtual void initData()
  {
    tools::FOR_EACH< detail::InitData >::loop(baseProblems);
  }

  virtual void finalizeData()
  {
    tools::FOR_EACH< detail::FinalizeData >::loop(baseProblems);
  }

  /**
   * At first the initData method is called for all baseProblems, then
   * the problems are filled with operators and coupling operators as well as
   * boundary conditions are added.
   **/
  virtual void initTimeInterface()
  {
    initData();
    fillOperators();
    fillBoundaryConditions();
    finalizeData();
  }


  virtual void fillCouplingOperators() {}
  virtual void fillCouplingBoundaryConditions() {}

  virtual void fillOperators()
  {
    tools::FOR_EACH< detail::FillOperators >::loop(baseProblems);
    fillCouplingOperators();
  }

  virtual void fillBoundaryConditions()
  {
    tools::FOR_EACH< detail::FillBoundaryConditions >::loop(baseProblems);
    fillCouplingBoundaryConditions();
  }

  /// get the j-th solution-vector of the i-th problem
  template<int i>
  DOFVector<double> *getSolution(int j)
  { FUNCNAME("CouplingBaseProblem::getSolution<i>(j)");
    BOOST_STATIC_ASSERT_MSG(0 <= i && i < _LENGTH_<BaseProblemsTupleType>::value , "********** ERROR: BaseProblem-index out of range **********");

    TEST_EXIT(0 <= j && j <= _GET_<i>(baseProblems).getNumComponents())("Indices out of range!\n");
    return _GET_<i>(baseProblems).getSolution()->getDOFVector(j);
  }


  /// pointer to the j-th feSpace of the i-th problem
  template<int i>
  inline const FiniteElemSpace* getFeSpace(int j=0)
  { FUNCNAME("CouplingBaseProblem::getFeSpace<i>(j)");
    BOOST_STATIC_ASSERT_MSG(0 <= i && i < _LENGTH_<BaseProblemsTupleType>::value , "********** ERROR: BaseProblem index out of range **********");

    TEST_EXIT(0 <= j && j <= _GET_<i>(baseProblems).getNumComponents())("Indices out of range!\n");
    return _GET_<i>(baseProblems).getFeSpace(j);
  }

  std::string getName() const { return name; }


  ProblemType *getProblem(std::string name_)
  {
    ProblemType *prob = NULL;
    tools::FOR_EACH< detail::FindProblem >::loop1(baseProblems, name_, prob);
    if (prob)
      return prob;
    else
      throw(std::runtime_error("problem with given name '" + name_ + "' does not exist"));
  }

  template<typename BaseProblemType>
  BaseProblemType *getBaseProblem(std::string name_)
  {
    BaseProblemType *prob = NULL;
    tools::FOR_EACH< detail::FindBaseProblem >::loop1(baseProblems, name_, prob);
    if (prob)
      return prob;
    else
      throw(std::runtime_error("problem with given name '" + name_ + "' does not exist"));
  }

  // final overriders for some functions...

  virtual void serialize(std::ostream &out) {};
  virtual void deserialize(std::istream &in) {};

  // using CouplingIterationInterface::beginIteration
  virtual void beginIteration(AdaptInfo *adaptInfo) override
  {
    CouplingIterationInterface::beginIteration(adaptInfo);
  }

  // using CouplingIterationInterface::oneIteration
  virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo) override
  {
    return CouplingIterationInterface::oneIteration(adaptInfo,toDo);
  }

  // using CouplingIterationInterface::endIteration
  virtual void endIteration(AdaptInfo *adaptInfo) override
  {
    CouplingIterationInterface::endIteration(adaptInfo);
  }

  // using CProblemStat::getNumProblems
  virtual int getNumProblems() const override
  {
    return CProblemStat::getNumProblems();
  }

  // using CProblemStat::getProblem
  virtual ProblemType *getProblem(int number = 0) override
  {
    return CProblemStat::getProblem(number);
  }

protected:
  BaseProblemsTupleType baseProblems;

  unsigned dim;	// dimension of the meshes
  unsigned dow;	// dimension of the world

  std::string name;
};

} // end namespace extensions

} // end namespace AMDiS

#endif // COUPLING_BASE_PROBLEM_H
