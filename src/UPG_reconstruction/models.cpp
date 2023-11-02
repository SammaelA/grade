#include "models.h"
namespace upg
{

  void add_tri_data(upg::vec3 point, upg::vec3 n, upg::vec2 tex, UniversalGenMesh &mesh)
  {
    mesh.pos.push_back(point.x);
    mesh.pos.push_back(point.y);
    mesh.pos.push_back(point.z);
    mesh.norm.push_back(n.x);
    mesh.norm.push_back(n.y);
    mesh.norm.push_back(n.z);
    mesh.tc.push_back(tex.x);
    mesh.tc.push_back(tex.y);
  }
  void add_point_data(upg::vec3 point, UniversalGenMesh &mesh)
  {
    mesh.pos.push_back(point.x);
    mesh.pos.push_back(point.y);
    mesh.pos.push_back(point.z);
  }
}