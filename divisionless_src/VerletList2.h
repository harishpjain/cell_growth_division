#pragma once

#include <vector>
#include <algorithm>

#include <mpi14/mpi14.hpp>

#include "io/ElementFileWriter.h"

namespace AMDiS
{

  template <class ForwardIt, class T, class Compare = std::less<>>
  ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, Compare comp = {})
  {
      first = std::lower_bound(first, last, value, comp);
      return first != last && !comp(value, *first) ? first : last;
  }

  // Strategy:
  // 1. Send center and radius to all processors
  // 2. Each processor builds its own neighbour-list 
  //
  // Implemented using MPI_Allgather
  class VerletList
  {
    static constexpr int tag_overlap_ = 101;
    static constexpr int tag_verlet_list_ = 102;
    static constexpr int tag_size_list_ = 103;

  public:
    VerletList(Mesh* mesh, mpi14::Communicator comm)
      : comm_(comm)
      , rank_(comm.rank())
      , size_(comm.size())
      , buffer_(3*size_)
    {
      Parameters::get(mesh->getName() + "->dimension", dim_);
      
      std::string output = "", postfix = "";
      Parameters::get("output", output);
      Parameters::get("postfix", postfix);
      
      std::string new_filename1 = output + "/debug/overlap" + "_" + postfix + "_p" + std::to_string(rank_) + "_";
      Parameters::set("overlap->output->filename", new_filename1);
      std::string new_filename2 = output + "/debug/verlet" + "_" + postfix + "_p" + std::to_string(rank_) + "_";
      Parameters::set("verlet->output->filename", new_filename2);
      
      writer1_.reset(new io::ElementFileWriter("overlap", mesh, elementVector1_));
      writer2_.reset(new io::ElementFileWriter("verlet", mesh, elementVector2_));
    }

    void compute(AdaptInfo* adaptInfo, DOFVector<double> const& signedDist, double dist)
    {
      FUNCNAME("VerletList::compute()");
      MSG("start (VerletList::compute) = %20.16e\n", mpi14::now());
      
      Timer t;
      center_ = centerOfCell(signedDist);
      MSG("time (computeCenter) = %e\n", t.elapsed());
            
      t.reset();
      radius_ = 0.0;
      auto overlap = computeOverlap(signedDist, dist, center_, radius_);
      MSG("time (computeOverlap) = %e\n", t.elapsed());
	  
#ifdef DEBUG_OUTPUT
      elementVector1_.clear();
      for (int macroIndex : overlap)
        elementVector1_[macroIndex] = 1.0;
      writer1_->writeFiles(adaptInfo, true, 0);
#endif
      
#ifdef DEBUG_PERFORMANCE
      MPI_Barrier(comm_);
#endif
      
      // build an interactioni graph, by approximating cells by circles with radius and center
      t.reset();
      auto graph = exchangeGraph(center_, radius_);
      MSG("time (exchangeGraph) = %e\n", t.elapsed());
            
      // send marked macro element to all neighbours in graph
      t.reset();
      auto overlaps = exchangeOverlap(overlap, graph);
      MSG("time (exchangeOverlap) = %e\n", t.elapsed());
      
      // build verlet-list from overlaps-list
      t.reset();
      for (auto& o : overlaps)
        std::sort(o.second.begin(), o.second.end());

      verletList_.clear();
      for (int m : overlap) {
        for (auto o_it = overlaps.begin(); o_it != overlaps.end(); ++o_it) {
          int r = o_it->first;
          auto const& overlap2 = o_it->second;
          
          auto it = binary_find(overlap2.begin(), overlap2.end(), m);
          if (it != overlap2.end())
            verletList_.push_back(std::make_pair(r,m));
        }
      }
      MSG("time (buildVerletList) = %e\n", t.elapsed());
      
#ifdef DEBUG_OUTPUT
      elementVector2_.clear();
      for (auto rank_macro : verletList_)
        elementVector2_[rank_macro.second] = rank_macro.first + 1;
      writer2_->writeFiles(adaptInfo, true, 0);
#endif
    }

    std::vector<std::pair<int,int>> const& get() const { return verletList_; }
	
	void computeNoOverlap(AdaptInfo* adaptInfo, DOFVector<double> const& signedDist, double dist)
    {
      FUNCNAME("VerletList::computeNoOverlap()");
      MSG("start (VerletList::computeNoOverlap) = %20.16e\n", mpi14::now());
      
      Timer t;
      //center_ = centerOfCell(signedDist);
      //MSG("time (computeCenter) = %e\n", t.elapsed());
            
      t.reset();
      //radius_ = 0.0;
      auto noOverlap = computeMacroIndices(signedDist/*, dist/*, center_ ,radius_*/);
      MSG("time (computeOverlap) = %e\n", t.elapsed());
      
# ifdef DEBUG_OUTPUT
      elementVector1_.clear();
      for (int macroIndex : noOverlap)
        elementVector1_[macroIndex] = 1.0;
      writer1_->writeFiles(adaptInfo, true, 0);
#endif
      
#ifdef DEBUG_PERFORMANCE
      MPI_Barrier(comm_);
#endif
      
      // build a list of neighbours rank
      t.reset();
      auto graph = allNeighboursGraph();
      MSG("time (allNeighboursGraph) = %e\n", t.elapsed());
            
      // send marked macro element to all neighbours in graph
      t.reset();
      auto noOverlaps = exchangeOverlap(noOverlap, graph);
      MSG("time (exchangeOverlap) = %e\n", t.elapsed());
      
      // build verlet-list from overlaps-list
      t.reset();
      for (auto& o : noOverlaps)
        std::sort(o.second.begin(), o.second.end());

      verletListNoOverLap_.clear();
      for (int m : noOverlap) {
        for (auto o_it = noOverlaps.begin(); o_it != noOverlaps.end(); ++o_it) {
          int r = o_it->first;
          auto const& noOverlap2 = o_it->second;
          
          auto it = binary_find(noOverlap2.begin(), noOverlap2.end(), m);
          if (it != noOverlap2.end())
            verletListNoOverLap_.push_back(std::make_pair(r,m));
        }
      }
      MSG("time (buildVerletList) = %e\n", t.elapsed());
      
#ifdef DEBUG_OUTPUT
      elementVector2_.clear();
      for (auto rank_macro : verletList_)
        elementVector2_[rank_macro.second] = rank_macro.first + 1;
      writer2_->writeFiles(adaptInfo, true, 0);
#endif
    }

    std::vector<std::pair<int,int>> const& getNoOverlap() const { return verletListNoOverLap_; }

  private:

    /// Returns vector of MacroElement indices overlapping with interface of signed-distance function
    std::vector<int> computeOverlap(DOFVector<double> const& signedDist, double dist, WorldVector<double> const& center, double& radius) const
    {
      FUNCNAME("VerletList::computeOverlap()");
      MSG("start (VerletList::computeOverlap) = %20.16e\n", mpi14::now());
      
      std::vector<int> overlap;

      FiniteElemSpace const* feSpace = signedDist.getFeSpace();
      BasisFunction const* basisFcts = feSpace->getBasisFcts();
      int numBasisFct = basisFcts->getNumber();

      std::vector<DegreeOfFreedom> localIndices(numBasisFct);

      radius = 0.0;
      for (int macroIndex = 0; macroIndex < feSpace->getMesh()->getNumberOfMacros(); ++macroIndex) {
        TraverseStack stack;
        ElInfo *elInfo = stack.traverseFirstOneMacro(feSpace->getMesh(), macroIndex, -1, Mesh::CALL_LEAF_EL);

        while (elInfo) {
          Element *el = elInfo->getElement();

          basisFcts->getLocalIndices(el, feSpace->getAdmin(), localIndices);

          bool inElement = false;
          for (int i = 0; i < numBasisFct; ++i) {
            //if (std::abs(signedDist[localIndices[i]]) < dist) {
            if (signedDist[localIndices[i]] < dist) {
              inElement = true;
              break;
            }
          }

          if (inElement) {
            overlap.push_back(macroIndex);
            
            MacroElement* macroEl = elInfo->getMacroElement();
            for (int i = 0; i < elInfo->getMesh()->getDim()+1; ++i)
              radius = std::max(radius, distance(macroEl->getCoord(i), center));
            
            break;
          }

          elInfo = stack.traverseNext(elInfo);
        }
      }

      return overlap;
    }
	
    /// Returns vector of all MacroElement indices 
    std::vector<int> computeMacroIndices(DOFVector<double> const& signedDist/*, double dist/*, WorldVector<double> const& center, double& radius*/) const
    {
      FUNCNAME("VerletList::computeOverlap()");
      MSG("start (VerletList::computeOverlap) = %20.16e\n", mpi14::now());
      
      std::vector<int> overlap;

      FiniteElemSpace const* feSpace = signedDist.getFeSpace();
      BasisFunction const* basisFcts = feSpace->getBasisFcts();
      //int numBasisFct = basisFcts->getNumber();

      //std::vector<DegreeOfFreedom> localIndices(numBasisFct);

      //radius = 0.0;
	  for (int macroIndex = 0; macroIndex < feSpace->getMesh()->getNumberOfMacros(); ++macroIndex) {
		overlap.push_back(macroIndex);}
      return overlap;
    }
    
    std::map<int, std::vector<int>> exchangeOverlap(std::vector<int>& myOverlap, std::vector<int> const& neighbours)
    {      
      FUNCNAME("VerletList::exchangeOverlap()");
      MSG("start (VerletList::exchangeOverlap) = %20.16e\n", mpi14::now());
      
      std::map<int, std::vector<int>> overlaps;
      
      // send and receive overlap vectors
      std::vector<mpi14::Request> requests;
      for (int r : neighbours) {
        auto req = comm_.isend(myOverlap, r, tag_overlap_);
        req.free();
        
        requests.push_back( comm_.irecv(overlaps[r], r, tag_overlap_) );
      }
      mpi14::wait_all(requests.begin(), requests.end());

      return overlaps;
    }
    
    
    WorldVector<double> centerOfCell(DOFVector<double> const& signedDist) const
    {
      FUNCNAME("VerletList::centerOfCell()");
      MSG("start (VerletList::centerOfCell) = %20.16e\n", mpi14::now());
      
      DOFConstIterator<double> vecIter(&signedDist, USED_DOFS);
      
      double minDist = 1.e5;
      DegreeOfFreedom dof = -1;
      for(vecIter.reset(); !vecIter.end(); ++vecIter) {
        if (*vecIter < minDist) {
          minDist = *vecIter;
          dof = vecIter.getDOFIndex();
        }
      }
      
      WorldVector<double> center;
      signedDist.getFeSpace()->getMesh()->getDofIndexCoords(&dof, signedDist.getFeSpace(), center);
      return center;
    }
    
    
    std::vector<int> exchangeGraph(WorldVector<double> const& center, double radius)
    {
      FUNCNAME("VerletList::exchangeGraph()");
      MSG("start (VerletList::exchangeGraph) = %20.16e\n", mpi14::now());
      
      std::array<double,3> pos {center[0], center[1], radius};
      std::cout << "DEBUG: Before MPI_Allgather" << std::endl;
      MPI_Allgather(pos.data(), 3, MPI_DOUBLE, buffer_.data(), 3, MPI_DOUBLE, comm_);
      std::cout << "DEBUG: Before MPI_Allgather" << std::endl;
      std::vector<int> neighbours; neighbours.reserve(8);
      for (int r = 0; r < size_; ++r) {
        if (r == rank_)
          continue;
        
        WorldVector<double> x; 
        x[0] = buffer_[3*r + 0];
        x[1] = buffer_[3*r + 1];
        double radius_r = buffer_[3*r + 2];
        
        if (distance(x, center) < radius + radius_r)
          neighbours.push_back(r);
      }
      
      return neighbours;
    }
	
	std::vector<int> allNeighboursGraph(){
		FUNCNAME("VerletList::allNeighboursGraph()");
		MSG("start (VerletList::allNeighboursGraph) = %20.16e\n", mpi14::now());
		
		std::vector<int> neighbours;
		neighbours.reserve(size_);
		for(int r = 0; r < size_; ++r){
			if(r==rank_)
				continue;
			neighbours.push_back(r);
		}
		return neighbours;
	}

    // periodic distance measure
    double distance(WorldVector<double> const& a, WorldVector<double> const& b) const
    {
      double d = 0.0;
      for (int i = 0; i < Global::getGeo(WORLD); ++i) {
        double di = std::abs(a[i] - b[i]);
        di-= dim_[i] * std::round(di/dim_[i]);
        d += di*di;
      }
      return std::sqrt(d);
    }
    
  public:
    
    WorldVector<double> getCenter() const { return center_; }
    double getRadius() const { return radius_; }
    
    
  private:

    mpi14::Communicator comm_;

    int rank_ = 0;
    int size_ = 1;

    WorldVector<double> dim_;
    WorldVector<double> center_;
    double radius_;
    
    std::vector<std::pair<int,int>> verletList_; // { (rank, macroIndex)... }
    std::vector<std::pair<int,int>> verletListNoOverLap_; // { (rank, macroIndex)... }
    
    std::map<int, double> elementVector1_;
    std::unique_ptr<io::ElementFileWriter> writer1_;
    
    std::map<int, double> elementVector2_;
    std::unique_ptr<io::ElementFileWriter> writer2_;
    
    
    std::vector<double> buffer_;
  };

} // end namespace AMDiS
