#pragma once

#include <string>
#include <deque>
#include <vector>

#include "AMDiS_fwd.h"
#include "DOFVector.h"
#include "io/FileWriterInterface.h"

namespace AMDiS
{
  /// \brief Create a contour grid of the zero-levelset given by a DOFVector
  /// and writes this grid to a VTU file.
  class ContourWriter
      : public FileWriterInterface
  {
    using Super = FileWriterInterface;
    
  public:
    
    /// Constructor. Takes a name to control init-file parameters 
    /// and a DOFVector representing the level-set function
    ContourWriter(std::string name, DOFVector<double>* vec);

    /// Extract contour from level-set function 'vec_' and write it to vtu/xyz/csv files
    virtual void writeFiles(AdaptInfo *adaptInfo, bool force,
			    int level = -1,
			    Flag traverseFlag = Mesh::CALL_LEAF_EL,
			    bool (*writeElem)(ElInfo*) = NULL) override;
    
    
  protected:
    
    /// Write .vtu file and (optional) .pvd file of the contour
    void writeVtuFile(AdaptInfo* adaptInfo, std::string const& fn,
                      std::vector< std::vector<int> > const& connectivity, 
                      std::vector<WorldVector<double>> const& vertices);
    
    /// Write points of the cfontour to .xyz file (one for each timestep)
    void writeXYZFile(AdaptInfo* adaptInfo, std::string const& fn, 
                      std::vector<WorldVector<double>> const& vertices);
    
    /// Write polygonal chains representing the contour as one line in a csv file.
    // NOTE: 2D only
    void writeCsvFile(AdaptInfo *adaptInfo, std::string const& fn, 
                      std::vector< std::vector<int> > const& connectivity, 
                      std::vector<WorldVector<double>> const& vertices);
    
    void writeBinaryFile(AdaptInfo *adaptInfo, std::string const& fn, 
                         std::vector< std::vector<int> > const& connectivity, 
                         std::vector<WorldVector<double>> const& vertices);
    
    
  private:
    
    // Extract the contour of a level-set function, by traversing all edges and taking those edges, 
    // with values of different sign on the vertices, as cutting edges. Calculate the position of the
    // interface cut on the edge and add this point to a list of vertices. All edge-cuts inside an 
    // element form a surface element.
    void contour(std::vector< std::vector<int> >& connectivity, 
                 std::vector<WorldVector<double>>& vertices, 
                 int level = -1, Flag traverseFlag = Mesh::CALL_LEAF_EL);
      
      
    // Calculate connected components of the surface-grid
    std::vector< std::deque<int> > groups(std::vector< std::vector<int> > const& connectivity, 
                                          int nVertices);
     
    // Helper function used in \ref groups
    std::pair<int,int> next_segment(int v_last, int c_last, 
                                    std::vector< std::vector<int> > const& connectivity, 
                                    std::vector< std::vector<int> > const& elements_of_vertex);
    
      
  private:
    
    // Name of the write in the initfile.
    std::string name_;
    
    // Level-set function
    DOFVector<double>* vec_;
    
    bool writeParaViewFormat_ = false;
    bool writeParaViewAnimation_ = false;
    bool writeCsvFormat_ = false;
    bool writeXYZFormat_ = false;
    bool writeBinaryFormat_ = false;
    
    int precision_ = 10;
    
    // vector that stores time-points and filenames
    std::vector<std::pair<double, std::string>> paraviewAnimationFrames_;
  };

} // end namespace AMDiS
