#include <cmath>
#include "tree_node.h"

u_g::vec3 norm(u_g::vec3 v1, u_g::vec3 v2)//right and up -> at screen
{
  return dgen::normalize_with_default(dgen::cross(v1, v2), u_g::vec3(1,0,0));
}

void add_rect(u_g::vec3 point, u_g::vec3 v1, u_g::vec3 v2, SimpleMeshData &mesh)
{
  add_tri(point, v1, v2, mesh);
  add_tri(point + v1 + v2, -v1, -v2, mesh);
}

void add_tri(u_g::vec3 point, u_g::vec3 v1, u_g::vec3 v2, SimpleMeshData &mesh)
{
  u_g::vec3 n = norm(v1, v2);
  u_g::vec3 p1 = point + v1, p2 = point + v2;
  mesh.add_tri_data(point, n, {0, 0});
  mesh.add_tri_data(p1, n, {0, 0});
  mesh.add_tri_data(p2, n, {0, 0});
}