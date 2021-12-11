#include "ContourWriter.h"

#include <string>
#include <map>
#include <deque>
#include <vector>
#include <fstream>

#include "AdaptInfo.h"
#include "BasisFunction.h"
#include "DOFVector.h"
#include "Element.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Initfile.h"
#include "Traverse.h"
#include "io/detail/VtkWriter.h"

#include "IndexMap.h"
#include "SplineCurve.h"

namespace AMDiS {
    
/// Constructor. Takes a name to control init-file parameters 
/// and a DOFVector representing the level-set function
ContourWriter::ContourWriter(std::string name, DOFVector<double>* vec)
  : name_(name)
  , vec_(vec)
{
  Super::readParameters(name);
  
  Parameters::get(name + "->ParaView format", writeParaViewFormat_);
  Parameters::get(name + "->ParaView animation", writeParaViewAnimation_);
  Parameters::get(name + "->XYZ format", writeXYZFormat_);  
  Parameters::get(name + "->CSV format", writeCsvFormat_);
  Parameters::get(name + "->Binary format", writeBinaryFormat_);
  Parameters::get(name + "->precision", precision_);
  
  if (writeCsvFormat_) { 
    std::ofstream file(filename + ".csv", std::ios_base::out);
    file << "#time,numChains,lenChain0,lenChain1,lenChain2,lenChain3,x0,y0,x1,y0,...,x31,y31\n";
  }
}


void ContourWriter::writeFiles(AdaptInfo *adaptInfo, bool force, int level,
                               Flag traverseFlag, bool (*writeElem)(ElInfo*) )
{
  if (!Super::doWriteTimestep(adaptInfo, force))
    return;
  
  std::vector<std::vector<int>> connectivity;
  std::vector<WorldVector<double>> vertices;
  contour(connectivity, vertices, level, traverseFlag);
  
  std::string fn;
  Super::getFilename(adaptInfo, fn);
  
  if (writeParaViewFormat_)
    writeVtuFile(adaptInfo, fn, connectivity, vertices);
  
  if (writeXYZFormat_)
    writeXYZFile(adaptInfo, fn, vertices);
  
  if (writeCsvFormat_)
    writeCsvFile(adaptInfo, fn, connectivity, vertices);
  
  if (writeBinaryFormat_)
    writeBinaryFile(adaptInfo, fn, connectivity, vertices);
}


void ContourWriter::writeVtuFile(AdaptInfo* adaptInfo, std::string const& fn,
                                 std::vector< std::vector<int> > const& connectivity, 
                                 std::vector<WorldVector<double>> const& vertices)
{ 
  int num_elements = connectivity.size();
  int num_vertices = vertices.size();
  
  std::map<int,int> cellType{
    {1,1}, // vertex
    {2,3}, // line
    {3,5}, // triangle
    {4,9}  // quad
  };

  { std::ofstream file(fn + ".vtu", std::ios_base::out);
    file << std::setprecision(precision_) << std::scientific;
    
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << num_vertices << "\" NumberOfCells=\"" <<  num_elements << "\">\n";
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto const& v : vertices) {
      for (int i = 0; i < Global::getGeo(WORLD); ++i)
        file << ' ' << v[i];
      for (int i = Global::getGeo(WORLD); i < 3; ++i)
        file << " 0.0";
      file << '\n';
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    file << "      <Cells>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";
    for (std::size_t i = 0; i < connectivity.size(); ++i)
      file << " " << (i + 1) * connectivity[i].size() << '\n';
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";
    for (auto const& element : connectivity)
      file << " " << cellType[ element.size() ] << '\n';
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";
    for (auto const& element : connectivity) {
      for (auto const& v : element)
        file << ' ' << v;
      file << '\n';
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
  }
  MSG("ParaView file written to %s\n", (fn + ".vtu").c_str());
  
  // write .pvd file
  if (writeParaViewAnimation_)
    io::VtkWriter::detail::updateAnimationFile(adaptInfo, fn + ".vtu", &paraviewAnimationFrames_, filename + ".pvd");
}


void ContourWriter::writeXYZFile(AdaptInfo* adaptInfo, std::string const& fn, 
                                 std::vector<WorldVector<double>> const& vertices)
{ 
  int num_vertices = vertices.size();
  
  std::ofstream file(fn + ".xyz", std::ios_base::out);
  file << std::setprecision(precision_) << std::scientific;
  file << num_vertices << '\n';
  file << "time: " << adaptInfo->getTime() << "\n";
  int nr = 0;
  for (auto const& v : vertices) {
    file << nr++;
    for (int i = 0; i < Global::getGeo(WORLD); ++i)
      file << ' ' << v[i];
    for (int i = Global::getGeo(WORLD); i < 3; ++i)
      file << " 0.0";
    file << '\n';
  }
  MSG("XYZ file written to %s\n", (fn + ".xyz").c_str());
}

void ContourWriter::writeCsvFile(AdaptInfo *adaptInfo, std::string const& fn, 
                                 std::vector< std::vector<int> > const& connectivity, 
                                 std::vector<WorldVector<double>> const& vertices)
{      
  if (vertices.empty())
    return;
  
  auto chains = groups(connectivity, vertices.size());
  
  int nTotalPoints = 32;
  Parameters::get(name_ + "->number of points", nTotalPoints);
  
  // calculate length of all chains
  std::vector<double> lengths;
  lengths.reserve(chains.size());
  double totalLength = 0.0;
  for (auto const& chain : chains) {
    assert( chain.size() > 1 );
    double length = 0;
    
    int v0 = chain.front();
    auto it = chain.begin();
    for (++it; it != chain.end(); ++it) {
      int v1 = *it;
      
      length += norm( vertices[v1] - vertices[v0] );
      v0 = v1;
    }
    
    lengths.push_back(length);
    totalLength += length;
  }
  
  std::vector<int> permutation(chains.size());
  std::iota(permutation.begin(), permutation.end(), 0);
  std::sort(permutation.begin(), permutation.end(), [&lengths](int i, int j) { return lengths[i] < lengths[j]; });
  
  // assign number of points per chain
  std::vector<int> nPoints(chains.size());
  int nAssignedPoints = 0;
  for (int i : permutation) {
    double l = lengths[i] / totalLength;
    
    // minimum 3 points per chain
    // maximum so that nTotalPoints not exceeded
    int p = std::max(3, std::min(nTotalPoints - nAssignedPoints, int(std::round(l * nTotalPoints))));
    nAssignedPoints += p;
    
    nPoints[i] = p;
  }
  
  // largest chain gets remaining points
  if (nAssignedPoints < nTotalPoints)
    nPoints[permutation.back()] += nTotalPoints - nAssignedPoints;
  
  // extract points of polygonal chains.
  std::vector< std::vector<WorldVector<double>> > chain_points(chains.size());
  for (int i = 0; i < chains.size(); ++i) {
    auto const& chain = chains[i];
    auto& points = chain_points[i];
    points.reserve(nPoints[i]);
    
    SplineCurve spline(index_map(chain.begin(), chain.end(), vertices));
    
    for (int j = 0; j < nPoints[i]; ++j) {
      double t = double(j)/(nPoints[i]-1);
      points.push_back(spline(t));
    }
  }
  
  
  // time,nChains,len(chain[0]),len(chain[1]),...,p[0]_x,p[0]_y,p[1]_x,p[1]_y,...
  { std::ofstream file(filename + ".csv", std::ios_base::app);
    file << std::setprecision(precision_) << std::scientific;
    file << adaptInfo->getTime() << ',';
    file << chain_points.size();
    for (auto const& chain : chain_points)
      file << ',' << chain.size();
    for (int i = chain_points.size(); i < 4; ++i)
      file << ",0";
    
    for (auto const& chain : chain_points)
      for (auto const& p : chain)
        file << ',' << p[0] << ',' << p[1];
      
    file << '\n';
  }
  MSG("CSV file written to %s\n", (filename + ".csv").c_str());
}


void ContourWriter::writeBinaryFile(AdaptInfo *adaptInfo, std::string const& /*fn*/, 
                                    std::vector< std::vector<int> > const& connectivity, 
                                    std::vector<WorldVector<double>> const& vertices)
{
  std::ofstream file(filename + ".bin", std::ios_base::app | std::ios_base::binary);
  
  double time = adaptInfo->getTime();
  file.write(reinterpret_cast<char*>(&time), sizeof(time));
  
  std::int32_t dow = Global::getGeo(WORLD);
  std::int32_t num_elements = connectivity.size();
  std::int32_t num_vertices = vertices.size();
  file.write(reinterpret_cast<char*>(&dow), sizeof(dow));
  file.write(reinterpret_cast<char*>(&num_elements), sizeof(num_elements));
  file.write(reinterpret_cast<char*>(&num_vertices), sizeof(num_vertices));
  
  for (auto const& element : connectivity) {
    std::int32_t s = element.size();
    file.write(reinterpret_cast<char*>(&s), sizeof(s));
    for (std::int32_t index : element)
      file.write(reinterpret_cast<char*>(&index), sizeof(index));
  }
  
  for (WorldVector<double> const& vertex : vertices)
    for (int i = 0; i < vertex.getSize(); ++i)
      file.write(reinterpret_cast<char const*>(&vertex[i]), sizeof(vertex[i]));
    
  MSG("Binary file written to %s\n", (filename + ".bin").c_str());
}


void ContourWriter::contour(std::vector< std::vector<int> >& connectivity, 
                            std::vector<WorldVector<double>>& vertices, 
                            int level, Flag traverseFlag)
{        
  using Vertices = std::map<std::pair<DegreeOfFreedom, DegreeOfFreedom>, 
                            std::pair<WorldVector<double>, int>>;
  Vertices vertices_tmp;
  int num_vertices = 0;
  
  FiniteElemSpace const* feSpace = vec_->getFeSpace();
  const BasisFunction *basisFcts = feSpace->getBasisFcts();
  int nBasisFcts = basisFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
        
  TraverseStack stack;

  int dim = feSpace->getMesh()->getDim();
  int nEdges = feSpace->getMesh()->getGeo(EDGE);
  ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(), level, traverseFlag | Mesh::FILL_COORDS);
  while (elInfo) {
    Element* el = elInfo->getElement();
    basisFcts->getLocalIndices(el, feSpace->getAdmin(), localIndices);
    
    std::vector<int> indices; indices.reserve(nEdges-1);
    
    // Check whether current element is cut by the zero level set.
    for (int e_i = 0; e_i < nEdges; ++e_i) {          
      int v_0 = el->getVertexOfEdge(e_i, 0);
      int v_1 = el->getVertexOfEdge(e_i, 1);
      
      // edge is cut by zero-level set
      double left = (*vec_)[localIndices[v_0]], right = (*vec_)[localIndices[v_1]];
      if ( (left<0.0 && right>=0.0) || (left>=0.0 && right<0.0) ) {
        double alpha = std::abs(left) / (std::abs(left) + std::abs(right));
        
        WorldVector<double> p = elInfo->getCoord(v_0);
        p += alpha * (elInfo->getCoord(v_1) - elInfo->getCoord(v_0));
        
        auto erg = vertices_tmp.insert(std::make_pair(el->getEdge(e_i), std::make_pair(p,num_vertices)));
        if (erg.second)
          ++num_vertices;
    
        indices.push_back(erg.first->second.second);
      }
    }
    
    if (!indices.empty())
      connectivity.push_back(indices);
    
    elInfo = stack.traverseNext(elInfo);
  }

  vertices.resize(vertices_tmp.size());
  for (auto const& v : vertices_tmp)
    vertices[v.second.second] = v.second.first;
}
  
  
std::vector< std::deque<int> > ContourWriter::groups(std::vector< std::vector<int> > const& connectivity, int nVertices)
{        
  int nElements = connectivity.size();
  
  std::vector< std::vector<int> > elements_of_vertex(nVertices); // = nodes
  for (int c_i = 0; c_i < nElements; ++c_i) {
    elements_of_vertex[connectivity[c_i][0]].push_back(c_i);
    elements_of_vertex[connectivity[c_i][1]].push_back(c_i);
  }
  
  std::vector< std::deque<int> > groups;
  for (int v_i = 0; v_i < nVertices; ++v_i) {
    auto& elements = elements_of_vertex[v_i];
    if (elements.empty())
      continue; // ignore single vertices
      
    std::deque<int> group; 
    group.push_back(v_i);
    
    int v_next = v_i;
    int c_next = elements.front();
    elements.pop_back();
    while (c_next >= 0) {
      std::tie(v_next, c_next) = next_segment(v_next, c_next, connectivity, elements_of_vertex);
      group.push_back(v_next);
      elements_of_vertex[v_next].clear();
    }
    
    v_next = v_i;
    c_next = elements.empty() ? -1 : elements.front();
    elements.pop_back();
    while (c_next >= 0) {
      std::tie(v_next, c_next) = next_segment(v_next, c_next, connectivity, elements_of_vertex);
      group.push_front(v_next);
      elements_of_vertex[v_next].clear();
    }
      
    groups.push_back(group);
  }
  
  return groups;
}
  
  
std::pair<int,int> ContourWriter::next_segment(int v_last, int c_last, 
                                               std::vector< std::vector<int> > const& connectivity, 
                                               std::vector< std::vector<int> > const& elements_of_vertex)
{
  auto const& vertices = connectivity[c_last];
  int v_i = vertices[0] == v_last ? vertices[1] : vertices[0];
  
  auto const& elements = elements_of_vertex[v_i];
  int c_i = elements.size() < 2 ? -1 : (elements[0] == c_last ? elements[1] : elements[0]);
              
  return std::make_pair(v_i,c_i);
}


} // end namespace AMDiS
