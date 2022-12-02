#pragma once

#include "AMDiS.h"

namespace AMDiS {

  template <class Term, class T, class Assigner>
  void assign(Term term, DOFVector<T>* result, int macroIndex, Assigner assign)
  {
    GenericOperatorTerm<Term> ot(term);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();

    Mesh* mesh = result->getFeSpace()->getMesh();

    std::vector<bool> visited(result->getUsedSize(), false);

    const FiniteElemSpace* resultFeSpace = result->getFeSpace();
    const BasisFunction *basisFcts = resultFeSpace->getBasisFcts();
    int nBasisFcts = basisFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    TraverseStack stack;

    Flag flags = Mesh::CALL_LEAF_EL/* | Mesh::FILL_BOUND | Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA*/;

    ElInfo* elInfo = stack.traverseFirstOneMacro(mesh, macroIndex, -1, flags);
    while (elInfo) {
      term.initElement(&ot, elInfo, NULL, NULL, basisFcts);
      basisFcts->getLocalIndices(elInfo->getElement(), resultFeSpace->getAdmin(), localIndices);

      for (int i = 0; i < nBasisFcts; i++) {
        if (!visited[localIndices[i]]) {
          assign(term(i), (*result)[localIndices[i]]);
          visited[localIndices[i]] = true;
        }
      }
      elInfo = stack.traverseNext(elInfo);
    }
  }

  template <class Term, class T, class Assigner>
  void assign(Term term, DOFVector<T>* result, std::vector<int> macroElements, Assigner assign)
  {
    GenericOperatorTerm<Term> ot(term);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();

    Mesh* mesh = result->getFeSpace()->getMesh();

    std::vector<bool> visited(result->getUsedSize(), false);

    const FiniteElemSpace* resultFeSpace = result->getFeSpace();
    const BasisFunction *basisFcts = resultFeSpace->getBasisFcts();
    int nBasisFcts = basisFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    TraverseStack stack;

    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;

    for (int macroIndex : macroElements) {
      ElInfo* elInfo = stack.traverseFirstOneMacro(mesh, macroIndex, -1, flags);
      while (elInfo) {
        term.initElement(&ot, elInfo, NULL, NULL, basisFcts);
        basisFcts->getLocalIndices(elInfo->getElement(), resultFeSpace->getAdmin(), localIndices);

        for (int i = 0; i < nBasisFcts; i++) {
          if (!visited[localIndices[i]]) {
            assign(term(i), (*result)[localIndices[i]]);
            visited[localIndices[i]] = true;
          }
        }
        elInfo = stack.traverseNext(elInfo);
      }
    }
  }


  inline void transfer(DOFVector<double> const& source, DOFVector<double>& target, int macroIndex)
  {
    const FiniteElemSpace* sourceFeSpace = source.getFeSpace();
    const FiniteElemSpace* feSpace = target.getFeSpace();
    const BasisFunction *basisFcts = feSpace->getBasisFcts();
    const BasisFunction *sourceBasisFcts = sourceFeSpace->getBasisFcts();

    int nBasisFcts = basisFcts->getNumber();
    int nSourceBasisFcts = sourceBasisFcts->getNumber();

    target.set(0.0);

    std::vector<DegreeOfFreedom> myLocalIndices(nBasisFcts);
    ElementVector sourceLocalCoeffs(nSourceBasisFcts);

    assert(feSpace->getMesh() == sourceFeSpace->getMesh());
	target = source; //newly added line
	/*
    DimVec<double> *coords = NULL;
    TraverseStack stack;
    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    ElInfo *elInfo = stack.traverseFirstOneMacro(feSpace->getMesh(), macroIndex, -1, flags);

    while (elInfo) {
      Element *el = elInfo->getElement();

      basisFcts->getLocalIndices(el, feSpace->getAdmin(), myLocalIndices);
      source.getLocalVector(el, sourceLocalCoeffs);

      for (int i = 0; i < nBasisFcts; i++) {
        if (target[myLocalIndices[i]] == 0.0) {
          coords = basisFcts->getCoords(i);
          target[myLocalIndices[i]] = sourceBasisFcts->evalUh(*coords, sourceLocalCoeffs);
        }
      }
      elInfo = stack.traverseNext(elInfo);
    }*/
  }


  template <class BinaryFunction>
  void transfer(DOFVector<double> const& vec1, DOFVector<double> const& vec2, DOFVector<double>& target, int macroIndex, BinaryFunction f)
  {
    const FiniteElemSpace* feSpace1 = vec1.getFeSpace();
    const FiniteElemSpace* feSpace2 = vec2.getFeSpace();
    const FiniteElemSpace* feSpace = target.getFeSpace();

    const BasisFunction *basisFcts1 = feSpace1->getBasisFcts();
    const BasisFunction *basisFcts2 = feSpace2->getBasisFcts();
    const BasisFunction *basisFcts = feSpace->getBasisFcts();

    int nBasisFcts1 = basisFcts1->getNumber();
    int nBasisFcts2 = basisFcts2->getNumber();
    int nBasisFcts = basisFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    ElementVector localCoeffs1(nBasisFcts1);
    ElementVector localCoeffs2(nBasisFcts2);

    assert(feSpace->getMesh() == feSpace1->getMesh());
    assert(feSpace1->getMesh() == feSpace2->getMesh());

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    ElInfo *elInfo = stack.traverseFirstOneMacro(feSpace->getMesh(), macroIndex, -1, flags);

    while (elInfo) {
      Element *el = elInfo->getElement();

      basisFcts->getLocalIndices(el, feSpace->getAdmin(), localIndices);
      vec1.getLocalVector(el, localCoeffs1);
      vec2.getLocalVector(el, localCoeffs2);

      for (int i = 0; i < nBasisFcts; i++) {
        if (target[localIndices[i]] == 0.0) {
          coords = basisFcts->getCoords(i);
          target[localIndices[i]] += f( basisFcts1->evalUh(*coords, localCoeffs1), basisFcts2->evalUh(*coords, localCoeffs2));
        }
      }
      elInfo = stack.traverseNext(elInfo);
    }
  }

  template <class BinaryFunction>
  double transfer_mm(DOFVector<double> const& vec1, DOFVector<double> const& vec2, DOFVector<double>& target, int macroIndex, BinaryFunction f)
  {
    const FiniteElemSpace* feSpace1 = vec1.getFeSpace();
    const FiniteElemSpace* feSpace2 = vec2.getFeSpace();
    const FiniteElemSpace* feSpace = target.getFeSpace();
    
    TEST_EXIT(feSpace1 == feSpace)("FeSpaces for vec1 and target must match!\n");

    const BasisFunction *basisFcts2 = feSpace2->getBasisFcts();
    const BasisFunction *basisFcts = feSpace->getBasisFcts();

    int nBasisFcts2 = basisFcts2->getNumber();
    int nBasisFcts = basisFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    ElementVector localCoeffs2(nBasisFcts2);
  
    DimVec<double> coords2(feSpace2->getMesh()->getDim(), NO_INIT);
    DualTraverse dualStack;
    ElInfo *elInfo1, *elInfo2;
    ElInfo *elInfoSmall, *elInfoLarge;
    WorldVector<double> worldVec;

    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse = dualStack.traverseFirstOneMacro(
      target.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags, flags,
      &elInfo1, &elInfo2, &elInfoSmall, &elInfoLarge);
    
    double interaction = 0.0;

    while (nextTraverse) {
      basisFcts->getLocalIndices(elInfo1->getElement(), feSpace->getAdmin(), localIndices);
      vec2.getLocalVector(elInfo2->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts; ++i) {
        if (target[localIndices[i]] == 0.0) {
          elInfo1->coordToWorld(*(basisFcts->getCoords(i)), worldVec);
          elInfo2->worldToCoord(worldVec, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target[localIndices[i]] += f(vec1[localIndices[i]], basisFcts2->evalUh(coords2, localCoeffs2));
            interaction +=   target[localIndices[i]];
          }  
        }
      }

      nextTraverse =
        dualStack.traverseNext(&elInfo1, &elInfo2, &elInfoSmall, &elInfoLarge);
    }
    return interaction;
  }
  
  template <class UnaryFunction>
  double transfer_mm(DOFVector<double> const& source, DOFVector<double>& target, int macroIndex, UnaryFunction f)
  {
    const FiniteElemSpace* sourcefeSpace = source.getFeSpace();
    const FiniteElemSpace* targetfeSpace = target.getFeSpace();
    
    const BasisFunction *sourceBasis = sourcefeSpace->getBasisFcts();
    const BasisFunction *targetBasis = targetfeSpace->getBasisFcts();

    int nBasisFctsSource = sourceBasis->getNumber();
    int nBasisFctsTarget = targetBasis->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasisFctsTarget);
    ElementVector localCoeffsSource(nBasisFctsSource);
  
    DimVec<double> coordsSource(sourcefeSpace->getMesh()->getDim(), NO_INIT);
    DualTraverse dualStack;
    ElInfo *elInfo1, *elInfo2;
    ElInfo *elInfoSmall, *elInfoLarge;
    WorldVector<double> worldVec;

    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse = dualStack.traverseFirstOneMacro(
      target.getFeSpace()->getMesh(), source.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags, flags,
      &elInfo1, &elInfo2, &elInfoSmall, &elInfoLarge);
    
    double interaction = 0.0;

    while (nextTraverse) {
      targetBasis->getLocalIndices(elInfo1->getElement(), targetfeSpace->getAdmin(), localIndices);
      source.getLocalVector(elInfo2->getElement(), localCoeffsSource);

      for (int i = 0; i < nBasisFctsTarget; ++i) {
        if (target[localIndices[i]] == 0.0) {
          elInfo1->coordToWorld(*(targetBasis->getCoords(i)), worldVec);
          elInfo2->worldToCoord(worldVec, &coordsSource);

          bool isPositive = true;
          for (int j = 0; j < coordsSource.getSize(); ++j) {
            if (coordsSource[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target[localIndices[i]] += f(sourceBasis->evalUh(coordsSource, localCoeffsSource));
            interaction +=   target[localIndices[i]];
          }  
        }
      }

      nextTraverse =
        dualStack.traverseNext(&elInfo1, &elInfo2, &elInfoSmall, &elInfoLarge);
    }
    return interaction;
  }
} // end namespace AMDiS
