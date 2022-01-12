#pragma once

#include <boost/serialization/vector.hpp>

#include "DOFSerializer.h"
#include "ExpressionAssigner.h"

#include <mpi14/mpi14.hpp>

namespace AMDiS
{

  template <class T>
  bool compare(T lhs, T rhs) { return lhs == rhs; }

  template <class T, class A>
  bool compare(std::vector<T,A> const& lhs, std::vector<T,A> const& rhs)
  {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), true,
      [](bool l, bool r) { return l && r; }, [](T const& l, T const& r) { return compare(l, r); });
  }

  template <class K, class T, class C, class A>
  bool compare(std::map<K,T,C,A> const& lhs, std::map<K,T,C,A> const& rhs)
  {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), true,
      [](bool l, bool r) { return l && r; }, [](auto const& l, auto const& r) { return compare(l.second, r.second); });
  }

  inline bool compare(double lhs, double rhs) { return std::abs(lhs - rhs) < 1.e-10; }

  class MacroData
  {
  public:
    MacroData() = default;

    MacroData(int macroIndex, DOFVector<double> const& vec)
    {
      MeshStructure structureCode;
      structureCode.init(vec.getFeSpace()->getMesh(), macroIndex);
      codeSize_ = structureCode.getNumElements();
      code_ = structureCode.getCode();

      DOFSerializer dofSerializer(vec.getFeSpace()->getAdmin());
      dofSerializer.gather(macroIndex, &vec, values_, false);
	    //std::cerr << macroIndex << ": " << values_[0] << "\n";
	  //std::cerr << values_[0] << "\n";// values_ are okay
    }

    void apply(int macroIndex, DOFVector<double>& vec, RefinementManager* refinementManager)
    {
      MeshStructure structureCode;
      structureCode.init(code_, codeSize_);

      structureCode.fitMeshToStructure(vec.getFeSpace()->getMesh(), refinementManager, false, macroIndex);
	    //std::cerr << values_[0] << "\n";
      DOFSerializer dofSerializer(vec.getFeSpace()->getAdmin());
      dofSerializer.scatter(macroIndex, values_, &vec, false); //To Dennis: What are these "values_". How is this class even getting these values
	}

  private:
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive& ar, const unsigned int /*version*/)
    {
      ar & codeSize_;
      ar & code_;
      ar & values_;
    }

    friend bool compare(MacroData const& lhs, MacroData const& rhs)
    {
      return compare(lhs.codeSize_, rhs.codeSize_) && compare(lhs.code_, rhs.code_) && compare(lhs.values_, rhs.values_);
    }

    int codeSize_ = 0;
    std::vector<uint64_t> code_;
    std::vector<double> values_;
  };


  // Communicate data on macro-element between processors
  class MacroCommunicator
  {
    static constexpr int tag_macrodata_ = 201;

  public:
    MacroCommunicator(mpi14::Communicator comm, std::vector<std::pair<int,int>> const& verletList)
      : comm_(comm)
      , verletList_(verletList){}
   // {}
  /*
    MacroCommunicator(mpi14::Communicator comm, std::vector<std::pair<int,int>> const& verletList , std::vector<std::pair<int,int>> const& verletListNoOverlap)
      : comm_(comm)
      , verletList_(verletList)
	  , verletListNoOverlap_(verletListNoOverlap)
    {}
  */
    // collect data from DOFVector and send to nighbouring ranks
    int scatter(DOFVector<double> const& vec)
    {
      FUNCNAME("MacroCommunicator::scatter()");
      MSG("start (MacroCommunicator::scatter) = %20.16e\n", mpi14::now());
      
      auto macrovector = init(vec);    //what does init does: Maps the DOF values to the macro mesh
      
      // send MacroData for all macros to rank
      for (auto const& rank_macro : macrovector) {                
        int r = rank_macro.first;
        auto req = comm_.isend(rank_macro.second, r, tag_macrodata_);
        req.free();
      }

      return macrovector.size();
    }
	
    void sendRanks(std::vector<std::vector<int>> out_ranks, int numRanks, int source){
      FUNCNAME("MacroCommunicator::sendRanks()");
      MSG("start (MacroCommunicator::sendRanks) = %20.16e\n", mpi14::now());

      for (int r = 0; r < numRanks; ++r){
        auto req = comm_.isend(out_ranks[source], r, tag_macrodata_);
        req.free();
      }
    }

    void recvRanks(std::vector<std::vector<int>>& in_ranks, int numRanks, int target){
      FUNCNAME("MacroCommunicator::recvRanks()");
      MSG("start (MacroCommunicator::recvRanks) = %20.16e\n", mpi14::now());
      std::vector<mpi14::Request> requests;
      for (int r = 0; r < numRanks; ++r)
      {
        requests.push_back(comm_.irecv(in_ranks[r], r, tag_macrodata_));
      }
      mpi14::wait_all(requests.begin(), requests.end());
    }

    void sendVolume(std::vector<std::vector<double>> out_volume, int numRanks, int source){
      FUNCNAME("MacroCommunicator::sendVolume()");
      MSG("start (MacroCommunicator::sendVolume) = %20.16e\n", mpi14::now());

      for (int r = 0; r < numRanks; ++r){
        auto req = comm_.isend(out_volume[source], r, tag_macrodata_);
        req.free();
      }
    }

    void recvVolume(std::vector<std::vector<double>>& in_volume, int numRanks, int target){
      FUNCNAME("MacroCommunicator::recvVolume()");
      MSG("start (MacroCommunicator::recvVolume) = %20.16e\n", mpi14::now());
      std::vector<mpi14::Request> requests;
      for (int r = 0; r < numRanks; ++r)
      {
        requests.push_back(comm_.irecv(in_volume[r], r, tag_macrodata_));
      }
      mpi14::wait_all(requests.begin(), requests.end());
    }



    // receive data from neighbouring ranks and apply to DOFVector
    template <class Operation>
    void gather(DOFVector<double>& vec, Operation apply)
    {
      FUNCNAME("MacroCommunicator::gather()");
      MSG("start (MacroCommunicator::gather) = %20.16e\n", mpi14::now());
      
      std::size_t numRanks = rankToMacros_.size();      
      std::vector<mpi14::Request> requests;
      
      // receive data objects
      std::vector<std::vector<MacroData>> macrovector(numRanks);
      for (std::size_t i = 0; i < numRanks; ++i) {
        int r = rankToMacros_[i].first;
        requests.push_back( comm_.irecv(macrovector[i], r, tag_macrodata_) );
      }
      mpi14::wait_all_apply(requests.begin(), requests.end(), [begin=requests.begin(),apply,&macrovector,&vec,this](auto it)
      {
        std::size_t i = std::distance(begin, it);
        
        for (std::size_t j = 0; j < macrovector[i].size(); ++j) {
          int macroIndex = this->rankToMacros_[i].second[j];
          
          this->coarseningManager_.globalCoarsen(vec.getFeSpace()->getMesh(), -20);
          macrovector[i][j].apply(macroIndex, vec, &this->refinementManager_);
          
          apply(macroIndex, vec);
        }
      });
      
    }
	

	
	
    // receive data from neighbouring ranks, apply to DOFVector and return contributions
    // expect 'apply' to have return argument in the form of contribution/(rank and macro)
    template <class Operation>
    void gather_and_deliver(DOFVector<double>& vec, std::vector<std::pair<std::size_t,double>>& contributions, Operation apply)
    {
      FUNCNAME("MacroCommunicator::gather()");
      MSG("start (MacroCommunicator::gather) = %20.16e\n", mpi14::now());
      
      std::size_t numRanks = rankToMacros_.size();      
      std::vector<mpi14::Request> requests;
      std::vector<size_t> requestRanks;
        
      // receive data objects
      std::vector<std::vector<MacroData>> macrovector(numRanks);
      for (std::size_t i = 0; i < numRanks; ++i) {
        int r = rankToMacros_[i].first;
        requests.push_back( comm_.irecv(macrovector[i], r, tag_macrodata_) );
        requestRanks.push_back(r);
      }

      mpi14::wait_all_apply(requests.begin(), requests.end(), [&requestRanks,&contributions,begin=requests.begin(),apply,&macrovector,&vec,this](auto it)
      {
        std::size_t i = std::distance(begin, it);

        double value = 0.0;
        for (std::size_t j = 0; j < macrovector[i].size(); ++j) {
          int macroIndex = this->rankToMacros_[i].second[j];
          this->coarseningManager_.globalCoarsen(vec.getFeSpace()->getMesh(), -20);
          macrovector[i][j].apply(macroIndex, vec, &this->refinementManager_);
          value += apply(macroIndex, vec);
        }
        contributions.push_back(std::make_pair(requestRanks[i],value));
      });
    }



    auto const& rankToMacros() const { return rankToMacros_; }
	  auto const& rankToMacros2() const { return rankToMacros2_; }
    protected:

    std::map<int, std::vector<MacroData>> init(DOFVector<double> const& vec) //make my own init
    {
      std::map<int, std::vector<int>> rm;
      for (auto const& rank_macro : verletList_)
        rm[rank_macro.first].push_back(rank_macro.second);      

      rankToMacros_.resize(rm.size());
      std::copy(rm.begin(), rm.end(), rankToMacros_.begin());     //copy rm to rankToMacros

      std::map<int, MacroData> macrodata; // macro-index => macro-data
      for (auto const& rank_macro : verletList_) {
        if (macrodata.count(rank_macro.second) == 0)
          macrodata.emplace( std::make_pair(rank_macro.second, MacroData{rank_macro.second, vec}) );
      }

      std::map<int, std::vector<MacroData>> macrovector; // rank => {macro-data_i}
      for (auto const& rank_macros : rankToMacros_) {
        macrovector[rank_macros.first].reserve(rank_macros.second.size());
        for (int m : rank_macros.second)
          macrovector[rank_macros.first].push_back(macrodata[m]);
      }
      return macrovector;
    }
	/*
	//returns macrodata for all points even without overlap
	  std::map<int, std::vector<MacroData>> initAll(DOFVector<double> const& vec) //make my own init
    {
      std::map<int, std::vector<int>> rm;
      for (auto const& rank_macro : verletListNoOverlap_)
        rm[rank_macro.first].push_back(rank_macro.second);      

      rankToMacros2_.resize(rm.size());
      std::copy(rm.begin(), rm.end(), rankToMacros2_.begin());     //copy rm to rankToMacros

      std::map<int, MacroData> macrodata; // macro-index => macro-data
      for (auto const& rank_macro : verletListNoOverlap_) {
        if (macrodata.count(rank_macro.second) == 0)
          macrodata.emplace( std::make_pair(rank_macro.second, MacroData{rank_macro.second, vec}) );
      }

      std::map<int, std::vector<MacroData>> macrovector; // rank => {macro-data_i}
      for (auto const& rank_macros : rankToMacros2_) {
        macrovector[rank_macros.first].reserve(rank_macros.second.size());
        for (int m : rank_macros.second)
          macrovector[rank_macros.first].push_back(macrodata[m]);
      }
      return macrovector;
    }*/

  protected:
    mpi14::Communicator comm_;

    std::vector< std::pair<int,int> > const& verletList_; // rank -> macro
	  //std::vector< std::pair<int,int> > const& verletListNoOverlap_; // rank -> macro
    std::vector< std::pair<int, std::vector<int>> > rankToMacros_;
    std::vector< std::pair<int, std::vector<int>> > rankToMacros2_;
	
    RefinementManager2d refinementManager_;
    CoarseningManager2d coarseningManager_;
  };

} // end namespace AMDiS
