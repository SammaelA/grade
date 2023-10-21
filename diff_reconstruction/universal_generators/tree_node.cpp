#include <cmath>
#include "tree_node.h"

u_g::vec3 norm(u_g::vec3 v1, u_g::vec3 v2)//right and up -> at screen
{
  u_g::vec3 n = {v1.y * v2.z - v1.z * v1.y, 
                 v1.z * v2.x - v1.x * v1.z, 
                 v1.x * v2.y - v1.y * v1.x};
  my_float len = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
  if (len == 0) 
  {
    n.x = 1;
  }
  else
  {
    n.x /= len;
    n.y /= len;
    n.z /= len;
  }
  return ;
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