#include <vector>
#include "p_v_structs.h"

struct SimpleMeshData
{
  static const int elem_size = 8;
  std::vector<my_float> data;
  void add_tri_data(u_g::vec3 point, u_g::vec3 n, u_g::vec2 tex)
  {
    data.push_back(point.x);
    data.push_back(point.y);
    data.push_back(point.z);
    data.push_back(n.x);
    data.push_back(n.y);
    data.push_back(n.z);
    data.push_back(tex.x);
    data.push_back(tex.y);
  }
};