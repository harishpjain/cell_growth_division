#pragma once

#include "PhaseFieldProblem.h"
#include "NavierStokesProblem.h"

#include "base_problems/CouplingBaseProblem2_cxx11.h"

namespace AMDiS { namespace extensions {

  class Coupling
      : public CouplingBaseProblem<ProblemStat, base_problems::PhaseFieldProblem, base_problems::NavierStokesProblem>
  {
    using Super = CouplingBaseProblem<ProblemStat, base_problems::PhaseFieldProblem, base_problems::NavierStokesProblem>;
    
  public:
    Coupling(base_problems::PhaseFieldProblem& vecProb, base_problems::NavierStokesProblem& nsProb)
      : Super("main", vecProb, nsProb)
      , vecProb_(vecProb)
      , nsProb_(nsProb)
    {
      Parameters::get("main->v0", v0_);
    }
    
    virtual void initData() override
    {
      Super::initData();
      
      phi_.reset(new DOFVector<double>(nsProb_.getFeSpace(0), "phi"));
      force_.reset(new DOFVector<WorldVector<double>>(nsProb_.getFeSpace(0), "mu_grad_phi"));
      
      for (int i = 0; i < vecProb_.getNumProblems(); ++i) {
        phis_.emplace_back(new DOFVector<double>(nsProb_.getFeSpace(0), "phi_i"));
        mus_.emplace_back(new DOFVector<double>(nsProb_.getFeSpace(0), "mu_i"));
        velocitiesX_.emplace_back(new DOFVector<double>(vecProb_.getFeSpace(i,0), "velX_i"));
        velocitiesY_.emplace_back(new DOFVector<double>(vecProb_.getFeSpace(i,0), "velY_i"));
      }
      velocities_.resize(vecProb_.getNumProblems());
      
      for (int i = 0; i < vecProb_.getNumProblems(); ++i) {
        velocities_[i][0] = velocitiesX_[i].get();
        velocities_[i][1] = velocitiesY_[i].get();
      }
      
      nsProb_.setPhi(phi_.get());
      
      fileWriter_.reset(new FileWriter("phase->space->output", getMesh(), phi_.get()));
    }
    
    virtual void solveInitialProblem(AdaptInfo* adaptInfo) override
    {
      vecProb_.solveInitialProblem(adaptInfo);
      
      getCoarseningManager()->globalCoarsen(getMesh(), -20);
      adapt_global_mesh();
      
      local_to_global(0, phis_);
      make_phi(phis_, *phi_);
      
      local_to_global(1, mus_);
      make_force(phis_, mus_, *force_);
      
      nsProb_.solveInitialProblem(adaptInfo);
    }
    
    virtual void transferInitialSolution(AdaptInfo *adaptInfo) override
    {
      Super::transferInitialSolution(adaptInfo);
      fileWriter_->writeFiles(adaptInfo, false);
    }
    
    virtual void initTimestep(AdaptInfo* adaptInfo) override
    {
      Super::initTimestep(adaptInfo);
      
      #pragma omp parallel for
      for (int i = 0; i < vecProb_.getNumProblems(); ++i) {
        velocitiesX_[i]->interpol(nsProb_.getProblem()->getSolution(0));
        velocitiesY_[i]->interpol(nsProb_.getProblem()->getSolution(1));
      }
    }
    
    virtual void closeTimestep(AdaptInfo* adaptInfo) override
    {
      vecProb_.closeTimestep(adaptInfo);
      
      adapt_global_mesh();
      
      local_to_global(0, phis_);
      make_phi(phis_, *phi_);
      
      local_to_global(1, mus_);
      make_force(phis_, mus_, *force_);
      
      // refine/coarsen global mesh along interface
      extensions::RefinementExpression(getMesh()).refine(function_(vecProb_.indicator("global", 0.95), 
        valueOf(*phi_)
      ));
      
      nsProb_.closeTimestep(adaptInfo);
      fileWriter_->writeFiles(adaptInfo, false);
    }
    
    virtual void fillCouplingOperators() override
    {
		
      for (int i = 0; i < vecProb_.getNumProblems(); ++i)
        fillCouplingOperatorsVec(vecProb_.getProblem(i), i);
      
      for (int i = 0; i < dow_; ++i) {
        // force term in NS-equation
        Operator* opInterface = new Operator(nsProb_.getFeSpace(i));
        addZOT(opInterface, v0_*componentOf(*force_,i));
        nsProb_.getProblem()->addVectorOperator(opInterface, i);
      }
    }
    
    void fillCouplingOperatorsVec(ProblemStat* prob, int i) 
    {
      Operator* opAdvect = new Operator(prob->getFeSpace(0));
      addFOT(opAdvect, -valueOf(velocities_[i]), GRD_PSI);
      
      prob->addMatrixOperator(opAdvect, 0,0);
    }
    
  protected:
    
    Mesh* getMesh() { return nsProb_.getMesh(0); }
    Mesh* getMesh(int i) { return vecProb_.getMesh(i); }
    
    RefinementManager* getRefinementManager() { return nsProb_.getProblem()->getRefinementManager(); }
    CoarseningManager* getCoarseningManager() { return nsProb_.getProblem()->getCoarseningManager(); }
    
    void adapt_global_mesh()
    {
      MeshStructure meshStructure{}; meshStructure.init(getMesh());
      for (int i = 0; i < vecProb_.getNumProblems(); ++i) {
        MeshStructure m{}; m.init(getMesh(i));
        meshStructure.merge(&m);
      }
      
      meshStructure.fitMeshToStructure(getMesh(), getRefinementManager(), false, -1, true);
    }
    
    void local_to_global(int comp, std::vector<std::unique_ptr<DOFVector<double>>>& globalSolution)
    {
      assert( globalSolution.size() == vecProb_.getNumProblems() );
      
      // copy local vector to global vector in parallel
      #pragma omp parallel for
      for (int i = 0; i < vecProb_.getNumProblems(); ++i)
        globalSolution[i]->interpol(vecProb_.getProblem(i)->getSolution(comp));
    }
    
    void make_phi(std::vector<std::unique_ptr<DOFVector<double>>> const& globalSolution, 
                  DOFVector<double>& phi) const
    {
      assert( globalSolution.size() == vecProb_.getNumProblems() );
      
      Max<double> f{};
      
      // merge individual global vectors
      phi.set(-1.0);
      for (int i = 0; i < vecProb_.getNumProblems(); ++i)
        transformDOF(globalSolution[i].get(), &phi, &phi, &f);
    }
    
    void make_force(std::vector<std::unique_ptr<DOFVector<double>>> const& phi, 
                    std::vector<std::unique_ptr<DOFVector<double>>> const& mu, 
                    DOFVector<WorldVector<double>>& force) const
    {
      assert( phi.size() == mu.size() && phi.size() == vecProb_.getNumProblems() );
      
      force << valueOf(*mu[0]) * gradientOf(*phi[0]);
      for (int i = 1; i < vecProb_.getNumProblems(); ++i)
        force << valueOf(force) + valueOf(*mu[i]) * gradientOf(*phi[i]);
    }
    
  private:
    base_problems::PhaseFieldProblem& vecProb_;
    base_problems::NavierStokesProblem& nsProb_;
    
    std::unique_ptr<DOFVector<double>> phi_;
    std::unique_ptr<DOFVector<WorldVector<double>>> force_;
    
    std::vector<std::unique_ptr<DOFVector<double>>> phis_;
    std::vector<std::unique_ptr<DOFVector<double>>> mus_;
    
    std::vector<std::unique_ptr<DOFVector<double>>> velocitiesX_, velocitiesY_;
    std::vector<WorldVector<DOFVector<double>*>> velocities_;
    
    std::unique_ptr<FileWriter> fileWriter_;
    
    double v0_ = 1.0;
    
    int dow_ = Global::getGeo(WORLD);
  };

} } // end namespaces
