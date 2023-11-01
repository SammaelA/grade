#include <cmath>
#include "tree_node.h"
#include "common_utils/template_vectors.h"
namespace upg
{
  upg::vec3 norm(upg::vec3 v1, upg::vec3 v2)//right and up -> at screen
  {
    return dgen::normalize_with_default(dgen::cross(v1, v2), upg::vec3(1,0,0));
  }

  upg::vec3 norm(upg::vec3 v)//right and up -> at screen
  {
    return dgen::normalize_with_default(v, upg::vec3(1,0,0));
  }

  upg::mat43 get_any_rot_mat(upg::vec3 axis, my_float angle)
  {
    vec3 ax = norm(axis);
    my_float c = cos(angle), s = sin(angle), x = ax.x, y = ax.y, z = ax.z;
    vec3 e1 = {c + (1 - c) * x * x, (1 - c) * x * y - s * z, (1 - c) * x * z + s * y};
    vec3 e2 = {(1 - c) * x * y - s * z, c + (1 - c) * y * y, (1 - c) * y * z + s * x};
    vec3 e3 = {(1 - c) * x * z + s * y, (1 - c) * x * y - s * z, c + (1 - c) * z * z};
    vec3 t = {0, 0, 0};
    return get_mat43(e1, e2, e3, t);
  }

  void add_rect(upg::vec3 point, upg::vec3 v1, upg::vec3 v2, UniversalGenMesh &mesh)
  {
    add_tri(point, v1, v2, mesh);
    add_tri(point + v1 + v2, -v1, -v2, mesh);
  }

  void add_tri(upg::vec3 point, upg::vec3 v1, upg::vec3 v2, UniversalGenMesh &mesh)
  {
    //upg::vec3 n = norm(v1, v2);
    upg::vec3 p1 = point + v1, p2 = point + v2;
    /*mesh.add_tri_data(point, n, {0, 0});
    mesh.add_tri_data(p1, n, {0, 0});
    mesh.add_tri_data(p2, n, {0, 0});*/
    add_point_data(point, mesh);
    add_point_data(p1, mesh);
    add_point_data(p2, mesh);
  }
}