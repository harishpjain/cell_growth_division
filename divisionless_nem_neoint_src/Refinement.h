#pragma once

#include <string>

namespace AMDiS
{

  // refinement indicator based on phase-field function
  inline auto indicator(std::string prefix, double threshold = 0.95)
  {
    int interface_ref = 15;
    int bulk_ref = 10;
    int outer_ref = 5;
    
    Parameters::get(prefix + "->interface refinements", interface_ref);
    Parameters::get(prefix + "->bulk refinements", bulk_ref);
    Parameters::get(prefix + "->outer refinements", outer_ref);
    
    return [interface_ref,bulk_ref,outer_ref,threshold](double const& phi) -> int
    {
      return (phi > -threshold && phi < threshold) ? interface_ref : (phi >= threshold) ? bulk_ref : outer_ref;
    };
  }
  
  // refinement indicator based on signed-dist function
  inline auto indicator2(std::string prefix, double dist)
  {
    int interface_ref = 15;
    int bulk_ref = 10;
    int outer_ref = 5;
    
    Parameters::get(prefix + "->interface refinements", interface_ref);
    Parameters::get(prefix + "->bulk refinements", bulk_ref);
    Parameters::get(prefix + "->outer refinements", outer_ref);
    
    return [dist,interface_ref,bulk_ref,outer_ref](double const& d) -> int
    {
      int threshold = 10*dist;
      int diff = interface_ref - outer_ref;
      
      int ref0 = outer_ref;
      int ref1 = ref0 + int(diff * (threshold - d) / (threshold - dist));
      
      return (d < -dist) ? bulk_ref : (d > threshold) ? ref0 : (d < dist) ? interface_ref : ref1;
    };
  }
  inline auto indicator2b(std::string prefix, double dist, bool boost_refinement=false)
  {

    int interface_ref = 15;
    int bulk_ref = 10;
    int outer_ref = 5;
    
    Parameters::get(prefix + "->interface refinements", interface_ref);
    Parameters::get(prefix + "->bulk refinements", bulk_ref);
    Parameters::get(prefix + "->outer refinements", outer_ref);

    if(boost_refinement){
      bulk_ref = 10;
      interface_ref = 10;
      outer_ref = 0;
    }
    
    return [dist,interface_ref,bulk_ref,outer_ref](double const& d) -> int
    {
      int threshold = 10*dist;
      int diff = interface_ref - outer_ref;
      
      int ref0 = outer_ref;
      int ref1 = ref0 + int(diff * (threshold - d) / (threshold - dist));
      
      return (d < -dist) ? bulk_ref : (d > threshold) ? ref0 : (d < dist) ? interface_ref : ref1;
    };
  }  
  
    // refinement indicator based on phase-field function
  inline auto indicator3(std::string prefix, double threshold = 0.95)
  {
    int interface_ref = 12;
    int bulk_ref = 12;  
    int outer_ref = 3;
    
    //Parameters::get(prefix + "->interface refinements", interface_ref);
    //Parameters::get(prefix + "->bulk refinements", bulk_ref);
    Parameters::get(prefix + "->outer refinements", outer_ref);
    
    return [interface_ref,bulk_ref,outer_ref,threshold](double const& phi) -> int
    {
      return (phi > -threshold && phi < threshold) ? interface_ref : (phi >= threshold) ? bulk_ref : outer_ref;
    };
  }

    // refinement indicator based on signed-dist function
  inline auto indicator4(std::string prefix, double dist)
  {
    int interface_ref = 11;
    int bulk_ref = 11;
    int outer_ref = 0;
    
    //Parameters::get(prefix + "->interface refinements", interface_ref);
    //Parameters::get(prefix + "->bulk refinements", bulk_ref);
    //Parameters::get(prefix + "->outer refinements", outer_ref);
    
    return [dist,interface_ref,bulk_ref,outer_ref](double const& d) -> int
    {
      int threshold = 10*dist;
      int diff = interface_ref - outer_ref;
      
      int ref0 = outer_ref;
      int ref1 = ref0 + int(diff * (threshold - d) / (threshold - dist));
      
      return (d < -dist) ? bulk_ref : (d > threshold) ? ref0 : (d < dist) ? interface_ref : ref1;
    };
  }
  // calculation of grid-width
  inline double mesh_width(Mesh* mesh)
  {
    FixVec<WorldVector<double>, VERTEX> coords(mesh->getDim(), NO_INIT);

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    double minH = 1e15;
    
    while (elInfo) {
      coords = elInfo->getCoords();
      double h = 0.0;
      for (int i = 0; i < coords.getSize(); i++) {
        for (int j = 0; j < coords.getSize(); j++) {
	  if (i != j)
            h = std::max(h,norm(coords[i]-coords[j]));
        }
      }
      minH = std::min(h, minH);
      elInfo = stack.traverseNext(elInfo);
    }

    return minH;
  }


  inline void adaptTimestep(AdaptInfo* adaptInfo)
  {
    std::vector<double> timesteps;
    std::vector<int> barriers;
    Parameters::get("adapt->timestep sequence", timesteps);
    Parameters::get("adapt->number sequence", barriers);
    TEST_EXIT(timesteps.size() == barriers.size())("timesteps.size != barriers.size\n");

    std::partial_sum(barriers.begin(), barriers.end(), barriers.begin());

    for (std::size_t i = 0; i < barriers.size(); ++i) {
      if (adaptInfo->getTimestepNumber() < barriers[i]) {
        adaptInfo->setTimestep(timesteps[i]);
        break;
      }
    }
  }

} // end namespace AMDiS
