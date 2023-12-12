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

  
  template <class BinaryFunction1, class BinaryFunction2>
  double transfer_mm_dual_binary(DOFVector<double> const& vec1, DOFVector<double> const& vec2, DOFVector<double>& target_t1, DOFVector<double>& target_t2, int macroIndex, BinaryFunction1 f1, BinaryFunction2 f2)
  {
    const FiniteElemSpace* feSpace1 = vec1.getFeSpace();
    const FiniteElemSpace* feSpace2 = vec2.getFeSpace();
    const FiniteElemSpace* feSpace_t1 = target_t1.getFeSpace();
    const FiniteElemSpace* feSpace_t2 = target_t2.getFeSpace();
    
    TEST_EXIT(feSpace1 == feSpace_t1)("FeSpaces for vec1 and target must match!\n");

    const BasisFunction *basisFcts2 = feSpace2->getBasisFcts();
    const BasisFunction *basisFcts_t1 = feSpace_t1->getBasisFcts();
    const BasisFunction *basisFcts_t2 = feSpace_t2->getBasisFcts();

    int nBasisFcts2 = basisFcts2->getNumber();
    int nBasisFcts_t1 = basisFcts_t1->getNumber();
    int nBasisFcts_t2 = basisFcts_t2->getNumber();

    std::vector<DegreeOfFreedom> localIndices_t1(nBasisFcts_t1);
    std::vector<DegreeOfFreedom> localIndices_t2(nBasisFcts_t2);
    ElementVector localCoeffs2(nBasisFcts2);
  
    DimVec<double> coords2(feSpace2->getMesh()->getDim(), NO_INIT);

    DualTraverse dualStack_t1;
    ElInfo *elInfo1_t1, *elInfo2_t1;
    ElInfo *elInfoSmall_t1, *elInfoLarge_t1;
    WorldVector<double> worldVec_t1;

    Flag flags_t1 = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_t1 = dualStack_t1.traverseFirstOneMacro(
      target_t1.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags_t1, flags_t1,
      &elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);
    
    double interaction_t1 = 0.0;

    while (nextTraverse_t1) {
      basisFcts_t1->getLocalIndices(elInfo1_t1->getElement(), feSpace_t1->getAdmin(), localIndices_t1);
      vec2.getLocalVector(elInfo2_t1->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts_t1; ++i) {
        if (target_t1[localIndices_t1[i]] == 0.0) {
          elInfo1_t1->coordToWorld(*(basisFcts_t1->getCoords(i)), worldVec_t1);
          elInfo2_t1->worldToCoord(worldVec_t1, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target_t1[localIndices_t1[i]] += f1(vec1[localIndices_t1[i]], basisFcts2->evalUh(coords2, localCoeffs2));
            interaction_t1 +=   target_t1[localIndices_t1[i]];
          }  
        }
      }

      nextTraverse_t1 =
        dualStack_t1.traverseNext(&elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);
    }

    DualTraverse dualStack_t2;
    ElInfo *elInfo1_t2, *elInfo2_t2;
    ElInfo *elInfoSmall_t2, *elInfoLarge_t2;
    WorldVector<double> worldVec_t2;

    Flag flags_t2 = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_t2 = dualStack_t2.traverseFirstOneMacro(
      target_t2.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags_t2, flags_t2,
      &elInfo1_t2, &elInfo2_t2, &elInfoSmall_t2, &elInfoLarge_t2);
    
    double interaction_t2 = 0.0;

    while (nextTraverse_t2) {
      basisFcts_t2->getLocalIndices(elInfo1_t2->getElement(), feSpace_t2->getAdmin(), localIndices_t2);
      vec2.getLocalVector(elInfo2_t2->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts_t2; ++i) {
        if (target_t2[localIndices_t2[i]] == 0.0) {
          elInfo1_t2->coordToWorld(*(basisFcts_t2->getCoords(i)), worldVec_t2);
          elInfo2_t2->worldToCoord(worldVec_t2, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target_t2[localIndices_t2[i]] += f2(vec1[localIndices_t2[i]], basisFcts2->evalUh(coords2, localCoeffs2));
            interaction_t2 +=  target_t2[localIndices_t2[i]];
          }  
        }
      }

      nextTraverse_t2 =
        dualStack_t2.traverseNext(&elInfo1_t2, &elInfo2_t2, &elInfoSmall_t2, &elInfoLarge_t2);
    }
    
    return interaction_t1;
  }

  
  template <class BinaryFunction1, class BinaryFunction2>
  double transfer_mm_dual_binary_part2(DOFVector<double> const& vec1, DOFVector<double> const& vec2, DOFVector<double>& target_t1, DOFVector<double>& target_t2, int macroIndex, BinaryFunction1 f1, BinaryFunction2 f2)
  {
    const FiniteElemSpace* feSpace1 = vec1.getFeSpace();
    const FiniteElemSpace* feSpace2 = vec2.getFeSpace();
    const FiniteElemSpace* feSpace_t1 = target_t1.getFeSpace();
    
    TEST_EXIT(feSpace1 == feSpace_t1)("FeSpaces for vec1 and target must match!\n");

    const BasisFunction *basisFcts2 = feSpace2->getBasisFcts();
    const BasisFunction *basisFcts_t1 = feSpace_t1->getBasisFcts();

    int nBasisFcts2 = basisFcts2->getNumber();
    int nBasisFcts_t1 = basisFcts_t1->getNumber();

    std::vector<DegreeOfFreedom> localIndices_t1(nBasisFcts_t1);
    ElementVector localCoeffs2(nBasisFcts2);
  
    DimVec<double> coords2(feSpace2->getMesh()->getDim(), NO_INIT);

    DualTraverse dualStack_t1;
    ElInfo *elInfo1_t1, *elInfo2_t1;
    ElInfo *elInfoSmall_t1, *elInfoLarge_t1;
    WorldVector<double> worldVec_t1;

    Flag flags_t1 = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_t1 = dualStack_t1.traverseFirstOneMacro(
      target_t1.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags_t1, flags_t1,
      &elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);
    
    double interaction_t1 = 0.0;

    while (nextTraverse_t1) {
      basisFcts_t1->getLocalIndices(elInfo1_t1->getElement(), feSpace_t1->getAdmin(), localIndices_t1);
      vec2.getLocalVector(elInfo2_t1->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts_t1; ++i) {
        if (target_t1[localIndices_t1[i]] == 0.0) {
          elInfo1_t1->coordToWorld(*(basisFcts_t1->getCoords(i)), worldVec_t1);
          elInfo2_t1->worldToCoord(worldVec_t1, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target_t1[localIndices_t1[i]] += f1(vec1[localIndices_t1[i]], basisFcts2->evalUh(coords2, localCoeffs2));
            target_t2[localIndices_t1[i]] += f2(vec1[localIndices_t1[i]], basisFcts2->evalUh(coords2, localCoeffs2));
            interaction_t1 += target_t1[localIndices_t1[i]];
          }  
        }
      }

      nextTraverse_t1 =
        dualStack_t1.traverseNext(&elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);
    }

    return interaction_t1;
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




  
  template <class UnaryFunction1, class UnaryFunction2>
  double transfer_mm_dual_unary(DOFVector<double> const& vec1, DOFVector<double> const& vec2, DOFVector<double>& target1, DOFVector<double>& target2, int macroIndex, UnaryFunction1 f1, UnaryFunction2 f2)
  {
    const FiniteElemSpace* feSpace1 = vec1.getFeSpace();
    const FiniteElemSpace* feSpace2 = vec2.getFeSpace();
    const FiniteElemSpace* feSpace_t1 = target1.getFeSpace();
    const FiniteElemSpace* feSpace_t2 = target2.getFeSpace();
    
    TEST_EXIT(feSpace1 == feSpace_t1)("FeSpaces for vec1 and target1 must match!\n");
    TEST_EXIT(feSpace1 == feSpace_t2)("FeSpaces for vec1 and target2 must match!\n");

    const BasisFunction *basisFcts2 = feSpace2->getBasisFcts();
    const BasisFunction *basisFcts_t1 = feSpace_t1->getBasisFcts();
    const BasisFunction *basisFcts_t2 = feSpace_t2->getBasisFcts();

    int nBasisFcts2 = basisFcts2->getNumber();
    int nBasisFcts_t1 = basisFcts_t1->getNumber();
    int nBasisFcts_t2 = basisFcts_t2->getNumber();

    std::vector<DegreeOfFreedom> localIndices_t1(nBasisFcts_t1);
    std::vector<DegreeOfFreedom> localIndices_t2(nBasisFcts_t2);
    ElementVector localCoeffs2(nBasisFcts2);
  
    DimVec<double> coords2(feSpace2->getMesh()->getDim(), NO_INIT);
    DualTraverse dualStack_t1;
    ElInfo *elInfo1_t1, *elInfo2_t1;
    ElInfo *elInfoSmall_t1, *elInfoLarge_t1;
    WorldVector<double> worldVec_t1;

    DualTraverse dualStack_t2;
    ElInfo *elInfo1_t2, *elInfo2_t2;
    ElInfo *elInfoSmall_t2, *elInfoLarge_t2;
    WorldVector<double> worldVec_t2;

    Flag flags = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_t1 = dualStack_t1.traverseFirstOneMacro(
      target1.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags, flags,
      &elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);

    bool nextTraverse_t2 = dualStack_t2.traverseFirstOneMacro(
      target2.getFeSpace()->getMesh(), vec2.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags, flags,
      &elInfo1_t2, &elInfo2_t2, &elInfoSmall_t2, &elInfoLarge_t2);
    
    double interaction_1 = 0.0;
    double interaction_2 = 0.0;

    while (nextTraverse_t1) {
      basisFcts_t1->getLocalIndices(elInfo1_t1->getElement(), feSpace_t1->getAdmin(), localIndices_t1);
      vec2.getLocalVector(elInfo2_t1->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts_t1; ++i) {
        if (target1[localIndices_t1[i]] == 0.0) {
          elInfo1_t1->coordToWorld(*(basisFcts_t1->getCoords(i)), worldVec_t1);
          elInfo2_t1->worldToCoord(worldVec_t1, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target1[localIndices_t1[i]] += f1(basisFcts2->evalUh(coords2, localCoeffs2)); //note only takes phi_2
            interaction_1 +=   target1[localIndices_t1[i]];
          }  
        }
      }

      nextTraverse_t1 =
        dualStack_t1.traverseNext(&elInfo1_t1, &elInfo2_t1, &elInfoSmall_t1, &elInfoLarge_t1);
    }

    while (nextTraverse_t2) {
      basisFcts_t2->getLocalIndices(elInfo1_t2->getElement(), feSpace_t2->getAdmin(), localIndices_t2);
      vec2.getLocalVector(elInfo2_t2->getElement(), localCoeffs2);

      for (int i = 0; i < nBasisFcts_t2; ++i) {
        if (target2[localIndices_t2[i]] == 0.0) {
          elInfo1_t2->coordToWorld(*(basisFcts_t2->getCoords(i)), worldVec_t2);
          elInfo2_t2->worldToCoord(worldVec_t2, &coords2);

          bool isPositive = true;
          for (int j = 0; j < coords2.getSize(); ++j) {
            if (coords2[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target2[localIndices_t2[i]] += f2(basisFcts2->evalUh(coords2, localCoeffs2));  //note only takes phi_2
            interaction_2 +=   target2[localIndices_t2[i]];
          }  
        }
      }

      nextTraverse_t2 =
        dualStack_t2.traverseNext(&elInfo1_t2, &elInfo2_t2, &elInfoSmall_t2, &elInfoLarge_t2);
    }

    //return interaction_1, interaction_2;
    return 0;
  }
  
  
  template <class UnaryFunction1, class UnaryFunction2>
  double transfer_mm_dual_unary_part2(DOFVector<double> const& source, DOFVector<double>& target_1, DOFVector<double>& target_2, int macroIndex, UnaryFunction1 f_1, UnaryFunction2 f_2)
  {
    const FiniteElemSpace* sourcefeSpace = source.getFeSpace();
    const FiniteElemSpace* targetfeSpace_1 = target_1.getFeSpace();
    const FiniteElemSpace* targetfeSpace_2 = target_2.getFeSpace();//
    
    const BasisFunction *sourceBasis = sourcefeSpace->getBasisFcts();
    const BasisFunction *targetBasis_1 = targetfeSpace_1->getBasisFcts();
    const BasisFunction *targetBasis_2 = targetfeSpace_2->getBasisFcts();//

    int nBasisFctsSource = sourceBasis->getNumber();
    int nBasisFctsTarget_1 = targetBasis_1->getNumber();
    int nBasisFctsTarget_2 = targetBasis_2->getNumber();//

    std::vector<DegreeOfFreedom> localIndices_1(nBasisFctsTarget_1);
    std::vector<DegreeOfFreedom> localIndices_2(nBasisFctsTarget_2);//
    ElementVector localCoeffsSource(nBasisFctsSource);
  
    DimVec<double> coordsSource(sourcefeSpace->getMesh()->getDim(), NO_INIT);
    DualTraverse dualStack_1;
    ElInfo *elInfo1_1, *elInfo2_1;
    ElInfo *elInfoSmall_1, *elInfoLarge_1;
    WorldVector<double> worldVec_1;

    Flag flags_1 = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_1 = dualStack_1.traverseFirstOneMacro(
      target_1.getFeSpace()->getMesh(), source.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags_1, flags_1,
      &elInfo1_1, &elInfo2_1, &elInfoSmall_1, &elInfoLarge_1);

    DualTraverse dualStack_2;//
    ElInfo *elInfo1_2, *elInfo2_2;
    ElInfo *elInfoSmall_2, *elInfoLarge_2;
    WorldVector<double> worldVec_2;

    Flag flags_2 = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS;
    bool nextTraverse_2 = dualStack_2.traverseFirstOneMacro(
      target_2.getFeSpace()->getMesh(), source.getFeSpace()->getMesh(),
      macroIndex, -1, -1, flags_2, flags_2,
      &elInfo1_2, &elInfo2_2, &elInfoSmall_2, &elInfoLarge_2);
    
    double interaction_1 = 0.0;

    while (nextTraverse_1) {
      targetBasis_1->getLocalIndices(elInfo1_1->getElement(), targetfeSpace_1->getAdmin(), localIndices_1);
      source.getLocalVector(elInfo2_1->getElement(), localCoeffsSource);

      for (int i = 0; i < nBasisFctsTarget_1; ++i) {
        if (target_1[localIndices_1[i]] == 0.0) {
          elInfo1_1->coordToWorld(*(targetBasis_1->getCoords(i)), worldVec_1);
          elInfo2_1->worldToCoord(worldVec_1, &coordsSource);

          bool isPositive = true;
          for (int j = 0; j < coordsSource.getSize(); ++j) {
            if (coordsSource[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target_1[localIndices_1[i]] += f_1(sourceBasis->evalUh(coordsSource, localCoeffsSource));
            interaction_1 +=   target_1[localIndices_1[i]];
          }  
        }
      }

      nextTraverse_1 =
        dualStack_1.traverseNext(&elInfo1_1, &elInfo2_1, &elInfoSmall_1, &elInfoLarge_1);
    }

    double interaction_2 = 0.0;

    while (nextTraverse_2) {
      targetBasis_2->getLocalIndices(elInfo1_2->getElement(), targetfeSpace_2->getAdmin(), localIndices_2);
      source.getLocalVector(elInfo2_2->getElement(), localCoeffsSource);

      for (int i = 0; i < nBasisFctsTarget_2; ++i) {
        if (target_2[localIndices_2[i]] == 0.0) {
          elInfo1_2->coordToWorld(*(targetBasis_2->getCoords(i)), worldVec_2);
          elInfo2_2->worldToCoord(worldVec_2, &coordsSource);

          bool isPositive = true;
          for (int j = 0; j < coordsSource.getSize(); ++j) {
            if (coordsSource[j] < -1.0e-5) {
              isPositive = false;
              break;
            }
          }

          if (isPositive){
            target_2[localIndices_2[i]] += f_2(sourceBasis->evalUh(coordsSource, localCoeffsSource));
            interaction_2 +=   target_2[localIndices_2[i]];
          }  
        }
      }

      nextTraverse_2 =
        dualStack_2.traverseNext(&elInfo1_2, &elInfo2_2, &elInfoSmall_2, &elInfoLarge_2);
    }
    return 0;//return interaction;
  }
  
  
} // end namespace AMDiS
