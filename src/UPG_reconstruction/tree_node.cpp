#include <cmath>
#include "tree_node.h"
#include "generation.h"
#include "common_utils/template_vectors.h"
#include "autodiff/autodiff.h"
namespace upg
{ 
  //const int MESH_REPEATS = 2;

  //const int SPINS = 16;

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
    vec3 ax = normalize_unsafe(axis);
    my_float c = cos(angle), s = sin(angle), x = ax.x, y = ax.y, z = ax.z;
    vec3 e1 = {c + (1 - c) * x * x, (1 - c) * x * y + s * z, (1 - c) * x * z - s * y};
    vec3 e2 = {(1 - c) * x * y - s * z, c + (1 - c) * y * y, (1 - c) * y * z + s * x};
    vec3 e3 = {(1 - c) * x * z + s * y, (1 - c) * z * y - s * x, c + (1 - c) * z * z};
    vec3 t = {0, 0, 0};
    return get_mat43(e1, e2, e3, t);
  }

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
    unsigned child_cnt() const override
    {
      return 0;
    }
    bool add_child(GenNode *node) override { return false; }
    std::vector<const GenNode *> get_children() const override
    {
      std::vector<const GenNode *> ret;
      return ret;
    }
  };

  void FreeTriangleNode_apply(const my_float *in, my_float *out)
  {
    for (int i=0;i<9;i++)
      out[i] = in[i];
  }

  void TriPrismNode_apply(const my_float *in, my_float *out)
  {
    out[0] = in[0];
    out[1] = 0;
    out[2] = in[1];
    out[3] = in[2];
    out[4] = 0;
    out[5] = in[3];
    out[6] = in[4];
    out[7] = 0;
    out[8] = in[5];

    out[9 + 0] = in[0];
    out[9 + 1] = 1;
    out[9 + 2] = in[1];
    out[9 + 3] = in[2];
    out[9 + 4] = 1;
    out[9 + 5] = in[3];
    out[9 + 6] = in[4];
    out[9 + 7] = 1;
    out[9 + 8] = in[5];
    for (int i = 0; i < 3; ++i)
    {
      out[18 + 18 * i + 0] = in[(0 + 2 * i) % 6];
      out[18 + 18 * i + 1] = 0;
      out[18 + 18 * i + 2] = in[(1 + 2 * i) % 6];
      out[18 + 18 * i + 3] = in[(0 + 2 * i) % 6];
      out[18 + 18 * i + 4] = 1;
      out[18 + 18 * i + 5] = in[(1 + 2 * i) % 6];
      out[18 + 18 * i + 6] = in[(2 + 2 * i) % 6];
      out[18 + 18 * i + 7] = 0;
      out[18 + 18 * i + 8] = in[(3 + 2 * i) % 6];

      out[18 + 18 * i + 9 + 0] = in[(2 + 2 * i) % 6];
      out[18 + 18 * i + 9 + 1] = 1;
      out[18 + 18 * i + 9 + 2] = in[(3 + 2 * i) % 6];
      out[18 + 18 * i + 9 + 3] = in[(2 + 2 * i) % 6];
      out[18 + 18 * i + 9 + 4] = 0;
      out[18 + 18 * i + 9 + 5] = in[(3 + 2 * i) % 6];
      out[18 + 18 * i + 9 + 6] = in[(0 + 2 * i) % 6];
      out[18 + 18 * i + 9 + 7] = 1;
      out[18 + 18 * i + 9 + 8] = in[(1 + 2 * i) % 6];
    }
  }

  void RotateNode_apply(const my_float *in, my_float *out)
  {
    vec3 ax = {in[0], in[1], in[2]};
    mat43 matr = get_any_rot_mat(ax, in[3]);
    vec3 v = {in[4], in[5], in[6]};
    v = mulv(matr, v);
    out[0] = v.x;
    out[1] = v.y;
    out[2] = v.z;
  }

  void ComplexRotateNode_apply(const my_float *in, my_float *out)
  {
    for (int i=0;i<3;i++)
      out[i] = in[7 + i] - in[4 + i];
    my_float data[7] = {in[0], in[1], in[2], in[3], out[0], out[1], out[2]};
    RotateNode_apply(data, out);
    for (int i=0;i<3;i++)
      out[i] += in[4 + i];
  }

  void ComplexRotateRepeatNode_apply(const my_float *in, my_float *out)
  {
    for (int i=0;i<3;i++)
    {
      out[i] = in[6 + i];
    }
    for (int j = 1; j < MESH_REPEATS; ++j)
    {
      my_float data[10] = {in[0], in[1], in[2], 2.0 * PI * (my_float)j / (my_float)MESH_REPEATS, in[3], in[4], in[5], in[6], in[7], in[8]};
      ComplexRotateNode_apply(data, out + j * 3);
    }
  }

  void SpinNode_8_apply(const my_float *in, my_float *out)
  {
    my_float y = 1;
    for (int i = 1; i < 8; ++i)
    {
      int j = 0;
      for (my_float ang = 1; ang <= SPINS; ++ang)
      {
        my_float prev_angle = (ang - 1) * 2.0 * M_PI / my_float(SPINS);
        my_float angle = ang * 2.0 * M_PI / my_float(SPINS);
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 0] = cos(prev_angle) * in[i - 1];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 1] = (y - 1)/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 2] = sin(prev_angle) * in[i - 1];

        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 3] = cos(prev_angle) * in[i];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 4] = y/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 5] = sin(prev_angle) * in[i];

        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 6] = cos(angle) * in[i - 1];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 7] = (y - 1)/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 0 + 8] = sin(angle) * in[i - 1];

        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 0] = cos(angle) * in[i];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 1] = y/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 2] = sin(angle) * in[i];

        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 3] = cos(angle) * in[i - 1];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 4] = (y - 1)/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 5] = sin(angle) * in[i - 1];

        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 6] = cos(prev_angle) * in[i];
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 7] = y/8;
        out[(i - 1) * SPINS * 2 * 9 + j * 2 * 9 + 9 + 8] = sin(prev_angle) * in[i];

        ++j;
      }
      ++y;
    }
  }

  class FreeTriangleNode : public PrimitiveNode
  {
  public:
    FreeTriangleNode(unsigned id) : PrimitiveNode(id) { name = "FreeTriangle"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      UniversalGenMesh mesh;
      mesh.pos.resize(9);
      if (out_jac)
      {
        out_jac->resize(9,9);
        ENZYME_EVALUATE_WITH_DIFF(FreeTriangleNode_apply, 9, 9, p.data(), mesh.pos.data(), out_jac->data());
      }
      else
        FreeTriangleNode_apply(p.data(), mesh.pos.data());
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      for (int i = 0; i < 9; i++)
        params.push_back({0, -1.0f, 1.0f, ParameterType::DIFFERENTIABLE, "p_" + std::to_string(i)});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 9;
    }
  };

  class FigureNode : public PrimitiveNode
  {
  public:
    FigureNode(unsigned id) : PrimitiveNode(id) { name = "Figure"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      UniversalGenMesh mesh;
      //creating for example cube
      add_rect({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, mesh);
      add_rect({0, 0, 0}, {0, 0, 1}, {0, 1, 0}, mesh);
      add_rect({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, mesh);

      add_rect({1, 1, 1}, {-1, 0, 0}, {0, -1, 0}, mesh);
      add_rect({1, 1, 1}, {0, -1, 0}, {0, 0, -1}, mesh);
      add_rect({1, 1, 1}, {0, 0, -1}, {-1, 0, 0}, mesh);
      if (out_jac)
        out_jac->resize(mesh.pos.size(), 0);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() const override
    {
      return 0;
    }
  };

  class SpinNode_8 : public PrimitiveNode
  {
  public:
    SpinNode_8(unsigned id) : PrimitiveNode(id) { name = "Spin_8"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      my_float data[8];
      for (int i = 0; i < 8; ++i)
      {
        data[i] = p[i];
      }
      UniversalGenMesh mesh;
      int mesh_size = (8 - 1) * SPINS * 2 * 9;
      mesh.pos.resize(mesh_size);
      //creating spin mesh
      if (out_jac)
      {
        out_jac->resize(mesh_size, param_cnt());
        ENZYME_EVALUATE_WITH_DIFF(SpinNode_8_apply, param_cnt(), mesh_size, p.data(), mesh.pos.data(), out_jac->data());
      }
      else
        SpinNode_8_apply(p.data(), mesh.pos.data());
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      for (int i = 0; i < 8; i++)
        params.push_back({0, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "spline_" + std::to_string(i)});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 8;
    }
  };

  class SphereNode : public PrimitiveNode
  {
  public:
    SphereNode(unsigned id) : PrimitiveNode(id) { name = "Sphere"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      //creating for example sphere
      my_float data[8], j = 0;
      for (int i = 0; i < 8; ++i)
      {
        data[i] = sqrt(abs(49.0 / 256.0 - pow(j / 8.0 - 7.0 / 16.0, 2)));
        j += 1;
      }
      UniversalGenMesh mesh;
      int mesh_size = (8 - 1) * SPINS * 2 * 9;
      mesh.pos.resize(mesh_size);
      SpinNode_8_apply(data, mesh.pos.data());
      if (out_jac)
        out_jac->resize(mesh.pos.size(), 0);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() const override
    {
      return 0;
    }
  };

  class ConeNode : public PrimitiveNode
  {
  public:
    ConeNode(unsigned id) : PrimitiveNode(id) { name = "Cone"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      //creating for example cone
      my_float data[8], j = 7;
      for (int i = 0; i < 8; ++i)
      {
        data[i] = j / 16.0;
        j -= 1;
      }
      UniversalGenMesh mesh;
      int mesh_size = (8 - 1) * SPINS * 2 * 9;
      mesh.pos.resize(mesh_size);
      SpinNode_8_apply(data, mesh.pos.data());
      for (int ang = 1; ang <= SPINS; ++ang)
      {
        my_float prev_angle = (ang - 1) * 2.0 * M_PI / my_float(SPINS);
        my_float angle = ang * 2.0 * M_PI / my_float(SPINS);
        add_tri({0, 0, 0}, {cos(angle) * 7.0 / 16.0, 0, sin(angle) * 7.0 / 16.0}, {cos(prev_angle) * 7.0 / 16.0, 0, sin(prev_angle) * 7.0 / 16.0}, mesh);
      }
      if (out_jac)
        out_jac->resize(mesh.pos.size(), 0);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() const override
    {
      return 0;
    }
  };

  class CylinderNode : public PrimitiveNode
  {
  public:
    CylinderNode(unsigned id) : PrimitiveNode(id) { name = "Cylinder"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      //creating for example cylinder
      my_float data[8];
      for (int i = 0; i < 8; ++i)
      {
        data[i] = 0.5;
      }
      UniversalGenMesh mesh;
      int mesh_size = (8 - 1) * SPINS * 2 * 9;
      mesh.pos.resize(mesh_size);
      SpinNode_8_apply(data, mesh.pos.data());
      for (int ang = 1; ang <= SPINS; ++ang)
      {
        my_float prev_angle = (ang - 1) * 2.0 * M_PI / my_float(SPINS);
        my_float angle = ang * 2.0 * M_PI / my_float(SPINS);
        add_tri({0, 7.0 / 8.0, 0}, {cos(angle) * 0.5, 0, sin(angle) * 0.5}, {cos(prev_angle) * 0.5, 0, sin(prev_angle) * 0.5}, mesh);
      }
      for (int ang = 1; ang <= SPINS; ++ang)
      {
        my_float prev_angle = (ang - 1) * 2.0 * M_PI / my_float(SPINS);
        my_float angle = ang * 2.0 * M_PI / my_float(SPINS);
        add_tri({0, 0, 0}, {cos(angle) * 0.5, 0, sin(angle) * 0.5}, {cos(prev_angle) * 0.5, 0, sin(prev_angle) * 0.5}, mesh);
      }
      if (out_jac)
        out_jac->resize(mesh.pos.size(), 0);
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() const override
    {
      return 0;
    }
  };

  class TriPrismNode : public PrimitiveNode
  {
  public:
    TriPrismNode(unsigned id) : PrimitiveNode(id) { name = "TriPrism"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      //creating for example prism
      UniversalGenMesh mesh;
      mesh.pos.resize(9 * 8);
      if (out_jac)
      {
        out_jac->resize(mesh.pos.size(), 6);
        ENZYME_EVALUATE_WITH_DIFF(TriPrismNode_apply, 6, mesh.pos.size(), p.data(), mesh.pos.data(), out_jac->data());
      }
      else
      {
        TriPrismNode_apply(p.data(), mesh.pos.data());
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      for (int i = 0; i < 6; i++)
        params.push_back({0, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "triangle_point_" + std::to_string(i)});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 6;
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
    unsigned child_cnt() const override
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
    std::vector<const GenNode *> get_children() const override
    {
      std::vector<const GenNode *> ret;
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
    ScaleNode(unsigned id) : OneChildNode(id) { name = "Scale"; }

    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      UniversalGenJacobian child_jac;
      UniversalGenMesh mesh = child->apply(out_jac ? &child_jac : nullptr);
      if (out_jac)
      {
        /*
        G = [Sx 0  0
             0  Sy 0
             0  0  Sz]
        M = [G 0 ... 0
             0 G ... 0
             ...
             0 0 ... G]
        H(v) = [v.x 0  0
                0  v.y 0
                0  0  v.z]
        M1 = [H(v0)
              H(v1)
              ...
              H(vn)]
        M' = [M*J | M1]
        */
        out_jac->resize(mesh.pos.size(), child_jac.get_yn() + 3);
        for (int i=0;i<child_jac.get_yn();i++)
        {
          for (int j=0;j<child_jac.get_xn();j++)
          {
            out_jac->at(3+i, j) = p[j % 3]*child_jac.at(i,j);
          }
        } 
        for (int i=0;i < 3; i++)
        {
          for (int j=0;j<child_jac.get_xn();j++)
          {
            if (j % 3 == i)
              out_jac->at(i, j) = (j % 3 == i) ? mesh.pos[j] : 0; //scale back to 
          }
        }
      }
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
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_x"});
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_y"});
      params.push_back({1.0f, 0.1f, 10.0f, ParameterType::DIFFERENTIABLE, "scale_z"});
      return params;
    }
    unsigned param_cnt() const override
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
    MoveNode(unsigned id) : OneChildNode(id) { name = "Move"; }

    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      UniversalGenJacobian child_jac;
      UniversalGenMesh mesh = child->apply(out_jac ? &child_jac : nullptr);
            
      if (out_jac)
      {
        /*
        G = [1  0  0
             0  1  0
             0  0  1]
        M = [G 0 ... 0
             0 G ... 0
             ...
             0 0 ... G]
        H(v) = [1  0  0
                0  1  0
                0  0  1]
        M1 = [H(v0)
              H(v1)
              ...
              H(vn)]
        M' = [M*J | M1]
        */
        out_jac->resize(mesh.pos.size(), child_jac.get_yn() + 3);
        for (int i=0;i<child_jac.get_yn();i++)
          for (int j=0;j<child_jac.get_xn();j++)
            out_jac->at(3+i, j) = child_jac.at(i,j);
      
        for (int i=0;i < 3; i++)
        {
          for (int j=0;j<child_jac.get_xn();j++)
          {
            if (j % 3 == i)
              out_jac->at(i, j) = (j % 3 == i) ? 1 : 0;
          }
        }
      }
      //applying moving
      for (int i = 0; i + 3 <= mesh.pos.size(); i += 3)
      {
        mesh.pos[i] += p[MOVE_X];
        mesh.pos[i + 1] += p[MOVE_Y];
        mesh.pos[i + 2] += p[MOVE_Z];
      }
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_z"});
      return params;
    }
    unsigned param_cnt() const override
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
    RotateNode(unsigned id) : OneChildNode(id) { name = "Rotate"; }

    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      /*my_float axis_x = p.get();
      my_float axis_y = p.get();
      my_float axis_z = p.get();
      my_float angle = p.get();*/
      UniversalGenJacobian child_jac;
      UniversalGenMesh mesh = child->apply(out_jac ? &child_jac : nullptr);
      //applying rotating
      vec3 ax = {p[AX_X], p[AX_Y], p[AX_Z]};
      mat43 matr = get_any_rot_mat(ax, p[ANGLE]);
      if (out_jac)
      {
        out_jac->resize(child_jac.get_xn(),child_jac.get_yn() + 4);
        UniversalGenJacobian G;
        G.resize(child_jac.get_xn(),child_jac.get_xn());
        for (int i = 0; i < mesh.pos.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(mesh.pos[i]);
          x.push_back(mesh.pos[i + 1]);
          x.push_back(mesh.pos[i + 2]);
          UniversalGenJacobian tmp;
          tmp.resize(3, 7);
          ENZYME_EVALUATE_WITH_DIFF(RotateNode_apply, 7, 3, x.data(), mesh.pos.data() + i, tmp.data());
          for (int a = 0; a < 4; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              out_jac->at(a, i + b) = tmp.at(a, b);
              if (a < 3)
              {
                G.at(i + a, i + b) = tmp.at(4 + a, b);
              }
            }
          }
        }
        for (int i=0;i<child_jac.get_yn();i++)
        {
          for (int j=0;j<child_jac.get_xn();j++)
          {
            for (int u = 0; u < child_jac.get_xn();u++)
            {
              if (u == 0)
              {
                out_jac->at(4+i, j) = G.at(u, j)*child_jac.at(i,u);
              }
              else
              {
                out_jac->at(4+i, j) += G.at(u, j)*child_jac.at(i,u);
              }
            }
          }
        } 
      }
      else
      {
        for (int i = 0; i < mesh.pos.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(mesh.pos[i]);
          x.push_back(mesh.pos[i + 1]);
          x.push_back(mesh.pos[i + 2]);
          RotateNode_apply(x.data(), mesh.pos.data() + i);
        }
      }
      matr = dgen::transposedInverse3x3(matr);
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
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_x"});
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_y"});
      params.push_back({1.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_z"});
      params.push_back({0.0f, -2*PI, 2*PI, ParameterType::DIFFERENTIABLE, "angle"});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 4;
    }
  };

  class ComplexRotateNode : public OneChildNode
  {
    static constexpr int AX_X = 0;
    static constexpr int AX_Y = 1;
    static constexpr int AX_Z = 2;
    static constexpr int ANGLE = 3;
    static constexpr int OFF_X = 4;
    static constexpr int OFF_Y = 5;
    static constexpr int OFF_Z = 6;

  public:
    ComplexRotateNode(unsigned id) : OneChildNode(id) { name = "ComplexRotate"; }

    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      /*my_float axis_x = p.get();
      my_float axis_y = p.get();
      my_float axis_z = p.get();
      my_float angle = p.get();*/
      UniversalGenJacobian child_jac;
      UniversalGenMesh mesh = child->apply(out_jac ? &child_jac : nullptr);
      //applying rotating
      vec3 ax = {p[AX_X], p[AX_Y], p[AX_Z]};
      mat43 matr = get_any_rot_mat(ax, p[ANGLE]);
      if (out_jac)
      {
        out_jac->resize(child_jac.get_xn(),child_jac.get_yn() + 7);
        UniversalGenJacobian G;
        G.resize(child_jac.get_xn(),child_jac.get_xn());
        for (int i = 0; i < mesh.pos.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(mesh.pos[i]);
          x.push_back(mesh.pos[i + 1]);
          x.push_back(mesh.pos[i + 2]);
          UniversalGenJacobian tmp;
          tmp.resize(3, 10);
          ENZYME_EVALUATE_WITH_DIFF(ComplexRotateNode_apply, 10, 3, x.data(), mesh.pos.data() + i, tmp.data());
          for (int a = 0; a < 7; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              out_jac->at(a, i + b) = tmp.at(a, b);
              if (a < 3)
              {
                G.at(i + a, i + b) = tmp.at(7 + a, b);
              }
            }
          }
        }
        for (int i=0;i<child_jac.get_yn();i++)
        {
          for (int j=0;j<child_jac.get_xn();j++)
          {
            for (int u = 0; u < child_jac.get_xn();u++)
            {
              if (u == 0)
              {
                out_jac->at(7+i, j) = G.at(u, j)*child_jac.at(i,u);
              }
              else
              {
                out_jac->at(7+i, j) += G.at(u, j)*child_jac.at(i,u);
              }
            }
          }
        } 
      }
      else
      {
        for (int i = 0; i < mesh.pos.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(mesh.pos[i]);
          x.push_back(mesh.pos[i + 1]);
          x.push_back(mesh.pos[i + 2]);
          ComplexRotateNode_apply(x.data(), mesh.pos.data() + i);
        }
      }
      matr = dgen::transposedInverse3x3(matr);
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
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_x"});
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_y"});
      params.push_back({1.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_z"});
      params.push_back({0.0f, -2*PI, 2*PI, ParameterType::DIFFERENTIABLE, "angle"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_x"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_y"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_z"});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 7;
    }
  };

  class ComplexRotateRepeatNode : public OneChildNode
  {
    static constexpr int AX_X = 0;
    static constexpr int AX_Y = 1;
    static constexpr int AX_Z = 2;
    //static constexpr int ANGLE = 3;
    static constexpr int OFF_X = 3;
    static constexpr int OFF_Y = 4;
    static constexpr int OFF_Z = 5;

  public:
    ComplexRotateRepeatNode(unsigned id) : OneChildNode(id) { name = "ComplexRotateRepeat"; }

    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      /*my_float axis_x = p.get();
      my_float axis_y = p.get();
      my_float axis_z = p.get();
      my_float angle = p.get();*/
      UniversalGenJacobian child_jac;
      UniversalGenMesh mesh = child->apply(out_jac ? &child_jac : nullptr);
      //applying rotating
      vec3 ax = {p[AX_X], p[AX_Y], p[AX_Z]};
      mat43 matr[MESH_REPEATS];
      for (int q = 0; q < MESH_REPEATS; ++q)
      {
        matr[q] = get_any_rot_mat(ax, 2 * PI * q / MESH_REPEATS);
      }
      std::vector<my_float> old_mesh;
      old_mesh.insert(old_mesh.end(), mesh.pos.begin(), mesh.pos.end());
      mesh.pos.resize(mesh.pos.size() * MESH_REPEATS);
      if (out_jac)
      {
        out_jac->resize(child_jac.get_xn() * MESH_REPEATS,child_jac.get_yn() + 6);
        UniversalGenJacobian G;
        G.resize(child_jac.get_xn() * MESH_REPEATS,child_jac.get_xn());
        for (int i = 0; i < old_mesh.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(old_mesh[i]);
          x.push_back(old_mesh[i + 1]);
          x.push_back(old_mesh[i + 2]);
          UniversalGenJacobian tmp;
          tmp.resize(3 * MESH_REPEATS, 9);
          std::vector <my_float> data(3 * MESH_REPEATS);
          ENZYME_EVALUATE_WITH_DIFF(ComplexRotateRepeatNode_apply, 9, 3 * MESH_REPEATS, x.data(), data.data(), tmp.data());
          for (int a = 0; a < MESH_REPEATS; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              mesh.pos[a * old_mesh.size() + i + b] = data[a * 3 + b];
            }
          }
          for (int a = 0; a < 6; ++a)
          {
            for (int c = 0; c < MESH_REPEATS; ++c)
            {
              for (int b = 0; b < 3; ++b)
              {
                out_jac->at(a, i + c * old_mesh.size() + b) = tmp.at(a, c * 3 + b);
                if (a < 3)
                {
                  G.at(i + a, i + c * old_mesh.size() + b) = tmp.at(6 + a, c * 3 + b);
                }
              }
            }
          }
        }
        for (int i=0;i<child_jac.get_yn();i++)
        {
          for (int j=0;j<child_jac.get_xn() * MESH_REPEATS;j++)
          {
            for (int u = 0; u < child_jac.get_xn();u++)
            {
              if (u == 0)
              {
                out_jac->at(6+i, j) = G.at(u, j)*child_jac.at(i,u);
              }
              else
              {
                out_jac->at(6+i, j) += G.at(u, j)*child_jac.at(i,u);
              }
            }
          }
        } 
      }
      else
      {
        for (int i = 0; i < old_mesh.size(); i += 3)
        {
          std::vector<my_float> x;
          x.insert(x.end(), p.begin(), p.end());
          x.push_back(old_mesh[i]);
          x.push_back(old_mesh[i + 1]);
          x.push_back(old_mesh[i + 2]);

          std::vector <my_float> data(3 * MESH_REPEATS);
          ComplexRotateRepeatNode_apply(x.data(), data.data());
          for (int a = 0; a < MESH_REPEATS; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              mesh.pos[a * old_mesh.size() + i + b] = data[a * 3 + b];
            }
          }
        }
      }
      for (int q = 0; q < MESH_REPEATS; ++q)
      {
        matr[q] = dgen::transposedInverse3x3(matr[q]);
      }
      old_mesh.clear();
      old_mesh.insert(old_mesh.end(), mesh.norm.begin(), mesh.norm.end());
      mesh.norm.resize(mesh.norm.size() * MESH_REPEATS);
      for (int q = 0; q < MESH_REPEATS; ++q)
      {
        for (int i = 0; i + 3 <= old_mesh.size(); i += 3)
        {
          vec3 v = {old_mesh[i], old_mesh[i + 1], old_mesh[i + 2]};
          v = mulv(matr[q], v);
          mesh.norm[q * old_mesh.size() + i] = v.x;
          mesh.norm[q * old_mesh.size() + i + 1] = v.y;
          mesh.norm[q * old_mesh.size() + i + 2] = v.z;
        }
      }
      logerr("XXXX");
      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_x"});
      params.push_back({0.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_y"});
      params.push_back({1.0f, 0.0f, 1.0f, ParameterType::DIFFERENTIABLE, "rot_axis_z"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_x"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_y"});
      params.push_back({0.0f, -10.0f, 10.0f, ParameterType::DIFFERENTIABLE, "offset_z"});
      return params;
    }
    unsigned param_cnt() const override
    {
      return 6;
    }
  };

  class TwoChildNode : public TransformNode
  {
  protected:
    GenNode *left;
    GenNode *right;
  public:
    TwoChildNode(unsigned id) : TransformNode(id) { left = NULL; right = NULL; }
    unsigned child_cnt() const override
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
    std::vector<const GenNode *> get_children() const override
    {
      std::vector<const GenNode *> ret;
      ret.push_back(left);
      ret.push_back(right);
      return ret;
    }
  };

  class MergeNode : public TwoChildNode
  {
  public:
    MergeNode(unsigned id) : TwoChildNode(id) { name = "Merge"; }
    UniversalGenMesh  apply(UniversalGenJacobian *out_jac) override
    {
      UniversalGenJacobian ch1,ch2;
      UniversalGenMesh mesh1 = left->apply(out_jac ? &ch1 : nullptr);
      UniversalGenMesh mesh2 = right->apply(out_jac ? &ch2 : nullptr);
      //applying or
      UniversalGenMesh mesh = mesh1;
      mesh.pos.insert(mesh.pos.end(), mesh2.pos.begin(), mesh2.pos.end());
      mesh.norm.insert(mesh.norm.end(), mesh2.norm.begin(), mesh2.norm.end());
      mesh.tc.insert(mesh.tc.end(), mesh2.tc.begin(), mesh2.tc.end());

      if (out_jac)
      {
        out_jac->resize(ch1.get_xn() + ch2.get_xn(), ch1.get_yn() + ch2.get_yn());
        out_jac->clear();
        for (int i=0;i<ch1.get_yn();i++)
          for (int j=0;j<ch1.get_xn();j++)
            out_jac->at(i, j) = ch1.at(i,j);
        for (int i=0;i<ch2.get_yn();i++)
          for (int j=0;j<ch2.get_xn();j++)
            out_jac->at(ch1.get_yn() + i, ch1.get_xn() + j) = ch2.at(i,j);
      }

      return mesh;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      return params;
    }
    unsigned param_cnt() const override
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
        node = new FigureNode(id);
        break;
      case 6:
        node = new MergeNode(id);
        break;
      case 7:
        node = new SpinNode_8(id);
        break;
      case 8:
        node = new ComplexRotateNode(id);
        break;
      case 9:
        node = new ComplexRotateRepeatNode(id);
        break;
      case 10:
        node = new SphereNode(id);
        break;
      case 11:
        node = new ConeNode(id);
        break;
      case 12:
        node = new CylinderNode(id);
        break;
      case 13:
        node = new TriPrismNode(id);
        break;
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;

    }
    return node;
  }
}