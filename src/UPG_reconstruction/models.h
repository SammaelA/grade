#pragma once
#include <vector>
#include "p_v_structs.h"
namespace upg
{
  struct UniversalGenMesh
  {
    //triangle mesh pos.size()%9 == 0
    //norm and tc can be empty
    std::vector<float> pos; //vec3
    std::vector<float> norm; //vec3
    std::vector<float> tc; //vec2
  };

  void add_tri_data(upg::vec3 point, upg::vec3 n, upg::vec2 tex, UniversalGenMesh &mesh);
  void add_point_data(upg::vec3 point, UniversalGenMesh &mesh);
}