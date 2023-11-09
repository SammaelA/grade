#include <cmath>
#include "tree_node.h"
#include "common_utils/template_vectors.h"
namespace upg
{ 
  vec3 norm(upg::vec3 v1, upg::vec3 v2);
  vec3 norm(upg::vec3 v);
  void add_rect(upg::vec3 point, upg::vec3 v1, upg::vec3 v2, UniversalGenMesh &mesh);
  void add_tri(upg::vec3 point, upg::vec3 v1, upg::vec3 v2, UniversalGenMesh &mesh);
  mat43 get_any_rot_mat(upg::vec3 axis, my_float angle);

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
    vec3 e2 = {(1 - c) * x * y + s * z, c + (1 - c) * y * y, (1 - c) * y * z - s * x};
    vec3 e3 = {(1 - c) * x * z - s * y, (1 - c) * z * y + s * x, c + (1 - c) * z * z};
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

  class PrimitiveNode : public GenNode
  {
  public:
    PrimitiveNode(unsigned id) : GenNode(id) {}
    unsigned child_cnt() override
    {
      return 0;
    }
    bool add_child(GenNode *node) override { return false; }
    std::vector<GenNode *> childs() override
    {
      std::vector<GenNode *> ret;
      return ret;
    }
  };

  class FreeTriangleNode : public PrimitiveNode
  {
  public:
    FreeTriangleNode(unsigned id) : PrimitiveNode(id) { node_num = 0; name = "FreeTriangle"; }
    UniversalGenMesh apply() override
    {
      UniversalGenMesh mesh;
      for (int i=0;i<9;i++)
        mesh.pos.push_back(p[i]);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      for (int i = 0; i < 9; i++)
        params.push_back({0, -1.0f, 1.0f, ParameterType::DIFFERENTIABLE, "p_" + std::to_string(i)});
      return params;
    }
    unsigned param_cnt() override
    {
      return 9;
    }
  };

  class FigureNode : public PrimitiveNode
  {
  public:
    FigureNode(unsigned id) : PrimitiveNode(id) { node_num = 0; name = "Figure"; }
    UniversalGenMesh apply() override
    {
      UniversalGenMesh mesh;
      //creating for example cube
      add_rect({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, mesh);
      add_rect({0, 0, 0}, {0, 0, 1}, {0, 1, 0}, mesh);
      add_rect({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, mesh);

      add_rect({1, 1, 1}, {-1, 0, 0}, {0, -1, 0}, mesh);
      add_rect({1, 1, 1}, {0, -1, 0}, {0, 0, -1}, mesh);
      add_rect({1, 1, 1}, {0, 0, -1}, {-1, 0, 0}, mesh);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() override
    {
      return 0;
    }
  };

  template <int N>
  class SpinNode : public PrimitiveNode
  {
  public:
    SpinNode(unsigned id) : PrimitiveNode(id) { node_num = 1; name = "Spin"; }
    UniversalGenMesh apply() override
    {
      my_float data[N];
      for (int i = 0; i < N; ++i)
      {
        data[i] = p[i];
      }
      UniversalGenMesh mesh;
      //creating spin mesh
      my_float y = 1;
      for (int i = 1; i < N; ++i, ++y)
      {
        for (my_float ang = 1; ang <= 8; ++ang)
        {
          my_float prev_angle = (ang - 1) * M_PI / 4.0;
          my_float angle = ang * M_PI / 4.0;
          add_tri({cos(prev_angle) * data[i - 1], y - 1, sin(prev_angle) * data[i - 1]}, 
                   {cos(prev_angle) * data[i], y, sin(prev_angle) * data[i]},
                   {cos(angle) * data[i - 1], y - 1, sin(angle) * data[i - 1]}, mesh);
          add_tri({cos(angle) * data[i], y, sin(prev_angle) * data[i]}, 
                   {cos(angle) * data[i - 1], y - 1, sin(angle) * data[i - 1]},
                   {cos(prev_angle) * data[i], y, sin(prev_angle) * data[i]}, mesh);
        }
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      for (int i = 0; i < N; i++)
        params.push_back({0, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "spline_" + std::to_string(i)});
      return params;
    }
    unsigned param_cnt() override
    {
      return N;
    }
  };

  class TransformNode : public GenNode
  {
  public:
    TransformNode(unsigned id) : GenNode(id) {}
  };

  class OneChildNode : public TransformNode
  {
  protected:
    GenNode *child;
  public:
    OneChildNode(unsigned id) : TransformNode(id) { child = NULL; }
    unsigned child_cnt() override
    {
      return 1;
    }
    bool add_child(GenNode *node) override 
    {
      if (child == NULL)
      {
        child = node;
      }
      return false;
    }
    std::vector<GenNode *> childs() override
    {
      std::vector<GenNode *> ret;
      ret.push_back(child);
      return ret;
    }
  };

  class ScaleNode : public OneChildNode
  {
    static constexpr int SCALE_X = 0;
    static constexpr int SCALE_Y = 1;
    static constexpr int SCALE_Z = 2;

  public:
    ScaleNode(unsigned id) : OneChildNode(id) { node_num = 2; name = "Scale"; }

    UniversalGenMesh apply() override
    {
      /*my_float scale_x = p.get();
      my_float scale_y = p.get();
      my_float scale_z = p.get();*/

      UniversalGenMesh mesh = child->apply();
      //applying scaling
      for (int i = 0; i + 3 <= mesh.pos.size(); i += 3)
      {
        mesh.pos[i] *= p[SCALE_X];
        mesh.pos[i + 1] *= p[SCALE_Y];
        mesh.pos[i + 2] *= p[SCALE_Z];
        if (mesh.norm.size() >= i + 3 && (i + 3) % 9 == 0)
        {
          vec3 p1 = {mesh.pos[i - 6], mesh.pos[i - 5], mesh.pos[i - 4]};
          vec3 p2 = {mesh.pos[i - 3], mesh.pos[i - 2], mesh.pos[i - 1]};
          vec3 p3 = {mesh.pos[i], mesh.pos[i + 1], mesh.pos[i + 2]};
          vec3 n = upg::norm(p2 - p1, p3 - p1);
          for (int j = 0; j < 9; j += 3)
          {
            mesh.norm[i - 6 + j] = n.x;
            mesh.norm[i - 5 + j] = n.y;
            mesh.norm[i - 4 + j] = n.z;
          }
        }
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_x"});
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_y"});
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_z"});
      return params;
    }
    unsigned param_cnt() override
    {
      return 3;
    }
  };

  class MoveNode : public OneChildNode
  {
    static constexpr int MOVE_X = 0;
    static constexpr int MOVE_Y = 1;
    static constexpr int MOVE_Z = 2;

  public:
    MoveNode(unsigned id) : OneChildNode(id) { node_num = 3; name = "Move"; }

    UniversalGenMesh apply() override
    {
      //my_float move_x = p.diff_params[MOVE_X];
      /*my_float move_x = p.get();
      my_float move_y = p.get();
      my_float move_z = p.get();*/
      UniversalGenMesh mesh = child->apply();
      //applying moving
      for (int i = 0; i + 3 <= mesh.pos.size(); i += 3)
      {
        mesh.pos[i] += p[MOVE_X];
        mesh.pos[i + 1] += p[MOVE_Y];
        mesh.pos[i + 2] += p[MOVE_Z];
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({-50.0f, 0.0f, 50.0f, ParameterType::DIFFERENTIABLE, "move_x"});
      params.push_back({-50.0f, 0.0f, 50.0f, ParameterType::DIFFERENTIABLE, "move_y"});
      params.push_back({-50.0f, 0.0f, 50.0f, ParameterType::DIFFERENTIABLE, "move_z"});
      return params;
    }
    unsigned param_cnt() override
    {
      return 3;
    }
  };

  class RotateNode : public OneChildNode
  {
    static constexpr int AX_X = 0;
    static constexpr int AX_Y = 1;
    static constexpr int AX_Z = 2;
    static constexpr int ANGLE = 3;

  public:
    RotateNode(unsigned id) : OneChildNode(id) { node_num = 4; name = "Rotate"; }

    UniversalGenMesh apply() override
    {
      /*my_float axis_x = p.get();
      my_float axis_y = p.get();
      my_float axis_z = p.get();
      my_float angle = p.get();*/

      UniversalGenMesh mesh = child->apply();
      //applying rotating
      vec3 ax = {p[AX_X], p[AX_Y], p[AX_Z]};
      mat43 matr = get_any_rot_mat(ax, p[ANGLE]);
      for (int i = 0; i + 3 <= mesh.pos.size(); i += 3)
      {
        vec3 v = {mesh.pos[i], mesh.pos[i + 1], mesh.pos[i + 2]};
        v = mulv(matr, v);
        mesh.pos[i] = v.x;
        mesh.pos[i + 1] = v.y;
        mesh.pos[i + 2] = v.z;
      }
      for (int i = 0; i + 3 <= mesh.norm.size(); i += 3)
      {
        vec3 v = {mesh.norm[i], mesh.norm[i + 1], mesh.norm[i + 2]};
        v = mulv(matr, v);
        mesh.norm[i] = v.x;
        mesh.norm[i + 1] = v.y;
        mesh.norm[i + 2] = v.z;
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_x"});
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_y"});
      params.push_back({1.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_z"});
      params.push_back({0.0f, 0.0f, 6.284f, ParameterType::DIFFERENTIABLE, "angle"});
      return params;
    }
    unsigned param_cnt() override
    {
      return 4;
    }
  };

  class TwoChildNode : public TransformNode
  {
  protected:
    GenNode *left;
    GenNode *right;
  public:
    TwoChildNode(unsigned id) : TransformNode(id) { left = NULL; right = NULL; }
    unsigned child_cnt() override
    {
      return 2;
    }
    bool add_child(GenNode *node) override 
    {
      if (left == NULL)
      {
        left = node;
        return true;
      }
      else if (right == NULL)
      {
        right = node;
      }
      return false;
    }
    std::vector<GenNode *> childs() override
    {
      std::vector<GenNode *> ret;
      ret.push_back(left);
      ret.push_back(right);
      return ret;
    }
  };

  class AndNode : public TwoChildNode
  {
  public:
    AndNode(unsigned id) : TwoChildNode(id) { node_num = 5; name = "And"; }
    UniversalGenMesh apply() override
    {
      UniversalGenMesh mesh1 = left->apply();
      UniversalGenMesh mesh2 = right->apply();
      //applying and
      return mesh1;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() override
    {
      return 0;
    }
  };

  class OrNode : public TwoChildNode
  {
  public:
    OrNode(unsigned id) : TwoChildNode(id) { node_num = 6; name = "Or"; }
    UniversalGenMesh apply() override
    {
      UniversalGenMesh mesh1 = left->apply();
      UniversalGenMesh mesh2 = right->apply();
      //applying or
      UniversalGenMesh mesh = mesh1;
      mesh.pos.insert(mesh.pos.end(), mesh2.pos.begin(), mesh2.pos.end());
      mesh.norm.insert(mesh.norm.end(), mesh2.norm.begin(), mesh2.norm.end());
      mesh.tc.insert(mesh.tc.end(), mesh2.tc.begin(), mesh2.tc.end());
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() override
    {
      return 0;
    }
  };

  class SubtrNode : public TwoChildNode
  {
  public:
    SubtrNode(unsigned id) : TwoChildNode(id) { node_num = 7; name = "Subtract"; }
    UniversalGenMesh apply() override
    {
      UniversalGenMesh mesh1 = left->apply();
      UniversalGenMesh mesh2 = right->apply();
      //applying subtract
      return mesh1;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() override
    {
      return 0;
    }
  };
  
  GenNode *node_by_node_type_id(uint16_t num, unsigned id)
  {
    GenNode *node = NULL;
    switch(num)
    {
      case 1: 
        node = new FreeTriangleNode(id);
        break;
      case 2:
        node = new ScaleNode(id);
        break;
      case 3:
        node = new MoveNode(id);
        break;
      case 4:
        node = new RotateNode(id);
        break;
      case 5:
        node = new AndNode(id);
        break;
      case 6:
        node = new OrNode(id);
        break;
      case 7:
        node = new SubtrNode(id);
        break;
      default:
        node = new FigureNode(id);
        break;

    }
    return node;
  }
}