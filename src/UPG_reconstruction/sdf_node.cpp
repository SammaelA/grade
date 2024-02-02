#include "sdf_node.h"
#include <stdio.h>
#include "generation.h"
#include "autodiff/autodiff.h"

float __enzyme_autodiff(...);

namespace upg
{
  AABB SdfGenInstance::scene_bbox;
  class PrimitiveSdfNode : public SdfNode
  {
  public:
    PrimitiveSdfNode(unsigned id) : SdfNode(id) {}
    virtual bool add_child(SdfNode *node) override { return false; }
    virtual unsigned child_cnt() const override { return 0; }
    virtual std::vector<const SdfNode *> get_children() const override { return {}; }
  };

  static constexpr int MAX_PARAMS = 256;

  #define GET_DISTANCE_WITH_DIFF(func)                       \
  {                                                          \
    float p_args[MAX_PARAMS];                                \
    p_args[0] = pos.x;                                       \
    p_args[1] = pos.y;                                       \
    p_args[2] = pos.z;                                       \
    for (int i = 0; i < param_cnt(); i++)                    \
      p_args[3 + i] = p[i];                                  \
    if (ddist_dp)                                            \
    {                                                        \
      float d_p[MAX_PARAMS] = {0};                           \
      float d = __enzyme_autodiff((void *)func, p_args, d_p);\
      int offset = ddist_dp->size();                         \
      ddist_dp->resize(offset + param_cnt());                \
      for (int i = 0; i < param_cnt(); i++)                  \
        (*ddist_dp)[offset + i] = d_p[3 + i];                \
      (*ddist_dpos)[0] = d_p[0];                             \
      (*ddist_dpos)[1] = d_p[1];                             \
      (*ddist_dpos)[2] = d_p[2];                             \
    }                                                        \
    return func(p_args);                                     \
  }

  class OneChildSdfNode : public SdfNode
  {
  protected:
    SdfNode *child;
  public:
    OneChildSdfNode(unsigned id) : SdfNode(id) { child = NULL; }
    unsigned child_cnt() const override { return 1; }
    bool add_child(SdfNode *node) override 
    {
      if (child == NULL)
        child = node;
      return false;
    }
    std::vector<const SdfNode *> get_children() const override
    {
      return {child};
    }
  };

  class TwoChildSdfNode : public SdfNode
  {
  protected:
    SdfNode *left;
    SdfNode *right;
  public:
    TwoChildSdfNode(unsigned id) : SdfNode(id) { left = NULL; right = NULL; }
    unsigned child_cnt() const override
    {
      return 2;
    }
    bool add_child(SdfNode *node) override 
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
    std::vector<const SdfNode *> get_children() const override
    {
      return {left, right};
    }
  };


  inline float 
  diff_sphere_sdf(float params[4]) 
  {
    return std::max(1e-9f, sqrt(params[0] * params[0] + params[1] * params[1] + params[2] * params[2])) - params[3];
  }
    
  class SphereSdfNode : public PrimitiveSdfNode
  {
    static constexpr int RADIUS = 0;
  public:
    SphereSdfNode(unsigned id) : PrimitiveSdfNode(id) { name = "Sphere"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_sphere_sdf);
    }

    virtual unsigned param_cnt() const override { return 1; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      float max_r = 0.5*length(scene_bbox.max_pos-scene_bbox.min_pos);
      std::vector<ParametersDescription::Param> params;
      params.push_back({1,0.01f*max_r,max_r, ParameterType::DIFFERENTIABLE, "radius"});
      return params;
    }
    virtual AABB get_bbox() const override
    {
      float r = p[RADIUS];
      return AABB({-r,-r,-r}, {r,r,r});
    }
  };

  inline float
  diff_box_sdf(float params[6])
  {
    float q[3], max_q[3];
    q[0] = abs(params[0]) - params[3];
    q[1] = abs(params[1]) - params[4];
    q[2] = abs(params[2]) - params[5];

    max_q[0] = std::max(q[0], 0.f);
    max_q[1] = std::max(q[1], 0.f);
    max_q[2] = std::max(q[2], 0.f);

    float d1 = std::sqrt(std::pow(max_q[0], 2) + std::pow(max_q[1], 2) + std::pow(max_q[2], 2));
    float d2 = std::min(std::max(q[0],std::max(q[1], q[2])), 0.f);

    return d1 + d2;
  }

  class BoxSdNode : public PrimitiveSdfNode
  {
    static constexpr int SIZE_X = 0;
    static constexpr int SIZE_Y = 1;
    static constexpr int SIZE_Z = 2;

  public:
    BoxSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Box"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_box_sdf);
    }

    virtual unsigned param_cnt() const override { return 3; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      glm::vec3 size = scene_bbox.max_pos-scene_bbox.min_pos;
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "size_x"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "size_y"});
      params.push_back({5,0.01f*size.z,size.z, ParameterType::DIFFERENTIABLE, "size_z"});

      return params;
    }
    virtual AABB get_bbox() const override
    {
      return AABB({-p[SIZE_X], -p[SIZE_Y], -p[SIZE_Z]}, {p[SIZE_X], p[SIZE_Y], p[SIZE_Z]});
    }
  };

  inline float
  diff_round_box_sdf(float params[7])
  {
    float q[3], max_q[3];
    q[0] = abs(params[0]) - params[3];
    q[1] = abs(params[1]) - params[4];
    q[2] = abs(params[2]) - params[5];

    max_q[0] = std::max(q[0], 0.f);
    max_q[1] = std::max(q[1], 0.f);
    max_q[2] = std::max(q[2], 0.f);

    float d1 = std::sqrt(std::pow(max_q[0], 2) + std::pow(max_q[1], 2) + std::pow(max_q[2], 2));
    float d2 = std::min(std::max(q[0],std::max(q[1], q[2])), 0.f);

    return d1 + d2 - params[6];
  }

  class RoundBoxSdNode : public PrimitiveSdfNode
  {
    static constexpr int SIZE_X = 0;
    static constexpr int SIZE_Y = 1;
    static constexpr int SIZE_Z = 2;
    static constexpr int RADIUS = 3;

  public:
    RoundBoxSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "RoundBox"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_round_box_sdf);
    }

    virtual unsigned param_cnt() const override { return 4; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      glm::vec3 size = scene_bbox.max_pos-scene_bbox.min_pos;
      float max_size = std::max(size.x, std::max(size.y, size.z));
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "size_x"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "size_y"});
      params.push_back({5,0.01f*size.z,size.z, ParameterType::DIFFERENTIABLE, "size_z"});
      params.push_back({0.3,0.01f*max_size,0.2f*max_size, ParameterType::DIFFERENTIABLE, "radius"});

      return params;
    }
    virtual AABB get_bbox() const override
    {
      return AABB({-(p[SIZE_X]+p[RADIUS]), -(p[SIZE_Y]+p[RADIUS]), -(p[SIZE_Z]+p[RADIUS])},
                  { (p[SIZE_X]+p[RADIUS]),  (p[SIZE_Y]+p[RADIUS]),  (p[SIZE_Z]+p[RADIUS])});
    }
  };

  inline float
  diff_cylinder_sdf(float params[5])
  {
    float d[2];
    d[0] = sqrt(std::pow(params[0], 2) + std::pow(params[2], 2)) - params[4];
    d[1] = std::abs(params[1]) - params[3];

    float d1 = 0, d2 = 0;
    float max_d[2];

    max_d[0] = std::max(d[0], 0.f);
    max_d[1] = std::max(d[1], 0.f);

    d1 = std::min(std::max(d[0], d[1]), 0.f);
    d2 = sqrt(std::pow(max_d[0], 2) + std::pow(max_d[1], 2));

    return d1 + d2;
  }

  class CylinderSdNode : public PrimitiveSdfNode
  {
    static constexpr int HEIGHT = 0;
    static constexpr int RADIUS = 1;

  public:
    CylinderSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Cylinder"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_cylinder_sdf);
    }

    virtual unsigned param_cnt() const override { return 2; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});
      params.push_back({2,0.01,10, ParameterType::DIFFERENTIABLE, "radius"});

      return params;
    }

    virtual AABB get_bbox() const override
    {
      return AABB({-p[RADIUS], -p[HEIGHT], -p[RADIUS]}, {p[RADIUS], p[HEIGHT], p[RADIUS]});
    }
  };

  inline float
  diff_prism_sdf(float params[5])
  {
    float q[3];
    q[0] = abs(params[0]);
    q[1] = abs(params[1]);
    q[2] = abs(params[2]);

    float d = std::max(q[2] - params[4], std::max(q[0] * 0.866025f + params[1] * 0.5f, -params[1]) - params[3] * 0.5f);

    return d;
  }

  class Prism : public PrimitiveSdfNode
  {
    static constexpr int H1 = 0;
    static constexpr int H2 = 1;

  public:
    Prism(unsigned id) : PrimitiveSdfNode(id) { name = "Prism"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_prism_sdf);
    }

    virtual unsigned param_cnt() const override { return 2; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height1"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height2"});

      return params;
    }
    virtual AABB get_bbox() const override
    {
      float max_h = std::max(p[H1], p[H2]);
      //I have no idea what the exact formula look like
      return AABB({-max_h, -max_h, -max_h}, {max_h, max_h, max_h});
    }
  };

  inline float
  diff_cone_sdf(float params[6])
  {
    float q[2], w[2], a[2], b[2];
    
    q[0] = params[5] * params[3] / params[4];
    q[1] = -params[5];

    w[0] = sqrt(std::pow(params[0], 2) + std::pow(params[2], 2));
    w[1] = params[1];

    a[0] = w[0] - q[0] * glm::clamp((w[0] * q[0] + w[1] * q[1]) / (q[0] * q[0] + q[1] * q[1]), 0.f, 1.f);
    a[1] = w[1] - q[1] * glm::clamp((w[0] * q[0] + w[1] * q[1]) / (q[0] * q[0] + q[1] * q[1]), 0.f, 1.f);

    b[0] = w[0] - q[0] * glm::clamp(w[0] / q[0], 0.f, 1.f);
    b[1] = w[1] - q[1];

    float k = glm::sign(q[1]);
    float v = std::min(a[0] * a[0] + a[1] * a[1], b[0] * b[0] + b[1] * b[1]);
    float s = std::max(k * (w[0] * q[1] - w[1] * q[0]), k * (w[1] - q[1]));
    float d = sqrt(v) * glm::sign(s);

    return d;
  }

  class Cone : public PrimitiveSdfNode
  {
    static constexpr int C1 = 0;
    static constexpr int C2 = 1;
    static constexpr int HEIGHT = 2;

  public:
    Cone(unsigned id) : PrimitiveSdfNode(id) { name = "Cone"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      GET_DISTANCE_WITH_DIFF(diff_cone_sdf);
    }

    virtual unsigned param_cnt() const override { return 3; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      glm::vec3 size = scene_bbox.max_pos-scene_bbox.min_pos;
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "c1"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "c2"});
      params.push_back({5,0.01f*size.z,size.z, ParameterType::DIFFERENTIABLE, "height"});

      return params;
    }
    virtual AABB get_bbox() const override
    {
      float max_s = std::max(p[C1], sqrtf(1-SQR(p[C2])));
      float r = max_s*p[HEIGHT];
      //I have no idea what the exact formula look like
      return AABB({-r, -p[HEIGHT], -r}, {r, 0, r});
    }
  };

  class MoveSdfNode : public OneChildSdfNode
  {
    static constexpr int MOVE_X = 0;
    static constexpr int MOVE_Y = 1;
    static constexpr int MOVE_Z = 2;
  public:
    MoveSdfNode(unsigned id) : OneChildSdfNode(id) { name = "Move"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      int offset;
      if (ddist_dp)
      {
        //[<prev_params>,ddist/dmove_params,<child_params>]
        offset = ddist_dp->size();
        ddist_dp->resize(offset + 3);
      }
      float d = child->get_distance(pos - glm::vec3(p[MOVE_X], p[MOVE_Y], p[MOVE_Z]), ddist_dp, ddist_dpos);

      if (ddist_dp)
      {
        // f(p,x,y) = g(h(p,y),x)
        // df/dp = dg/d1 * dh/dp
        // df/dx = dg/dx
        // df/dy = dg/d1 * dh/dy
        // g = child->get_distance()
        // h(pos, move) = pos - move

        (*ddist_dp)[offset+0] = -(*ddist_dpos)[0];
        (*ddist_dp)[offset+1] = -(*ddist_dpos)[1];
        (*ddist_dp)[offset+2] = -(*ddist_dpos)[2];

        //(*ddist_dpos)[0] = (*ddist_dpos)[0];
        //(*ddist_dpos)[1] = (*ddist_dpos)[1];
        //(*ddist_dpos)[2] = (*ddist_dpos)[2];
      }

      return d;
    }
    virtual unsigned param_cnt() const override { return 3; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,scene_bbox.min_pos.x,scene_bbox.max_pos.x, ParameterType::DIFFERENTIABLE, "move_x"});
      params.push_back({0,scene_bbox.min_pos.y,scene_bbox.max_pos.y, ParameterType::DIFFERENTIABLE, "move_y"});
      params.push_back({0,scene_bbox.min_pos.z,scene_bbox.max_pos.z, ParameterType::DIFFERENTIABLE, "move_z"});
      return params;
    }
    virtual AABB get_bbox() const override
    {
      AABB ch_bbox = child->get_bbox();
      glm::vec3 sh = glm::vec3(p[MOVE_X], p[MOVE_Y], p[MOVE_Z]);
      return AABB(ch_bbox.min_pos+sh, ch_bbox.max_pos+sh);
    }
  };

  void get_any_rot_mat_to_vec(my_float x, my_float y, my_float z, float angle, my_float *out)
  {
    float c = cos(angle), s = sin(angle);
    glm::vec3 e1 = {c + (1 - c) * x * x, (1 - c) * x * y + s * z, (1 - c) * x * z - s * y};
    glm::vec3 e2 = {(1 - c) * x * y - s * z, c + (1 - c) * y * y, (1 - c) * y * z + s * x};
    glm::vec3 e3 = {(1 - c) * x * z + s * y, (1 - c) * z * y - s * x, c + (1 - c) * z * z};
    out[0] = e1.x * x + e1.y * y + e1.z * z;
    out[1] = e2.x * x + e2.y * y + e2.z * z;
    out[2] = e3.x * x + e3.y * y + e3.z * z;
  }

  void RotateSdfNode_apply(const my_float *in, my_float *out)
  {
    out[0] = in[3];
    out[1] = in[4];
    out[2] = in[5];
    get_any_rot_mat_to_vec(cos(in[0]) * cos(in[1]), sin(in[0]) * cos(in[1]), sin(in[1]), in[2], out);
  }

  class RotateSdfNode : public OneChildSdfNode
  {
    static constexpr int AXIS_ANG_XY = 0;
    static constexpr int AXIS_ANG_Z = 1;
    static constexpr int ANGLE = 2;
  public:
    RotateSdfNode(unsigned id) : OneChildSdfNode(id) { name = "Rotate"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      int offset;
      float d = 0;
      if (ddist_dp)
      {
        //[<prev_params>,ddist/dmove_params,<child_params>]
        offset = ddist_dp->size();
        ddist_dp->resize(offset + 4);
      }
      

      if (ddist_dp)
      {
        // f(p,x,y) = g(h(p,y),x)
        // df/dp = dg/d1 * dh/dp
        // df/dx = dg/dx
        // df/dy = dg/d1 * dh/dy
        // g = child->get_distance()
        // h(pos, move) = matr * pos

        std::vector<float> x;
        x.insert(x.end(), p.begin(), p.end());
        x.push_back(pos.x);
        x.push_back(pos.y);
        x.push_back(pos.z);
        upg::UniversalGenJacobian tmp;
        tmp.resize(3, 7);
        std::vector<float> y(3);
        ENZYME_EVALUATE_WITH_DIFF(RotateSdfNode_apply, 7, 3, x.data(), y.data(), tmp.data());
        glm::vec3 ps = {y[0], y[1], y[2]};
        d = child->get_distance(ps, ddist_dp, ddist_dpos);

        float X = (*ddist_dpos)[0], Y = (*ddist_dpos)[1], Z = (*ddist_dpos)[2];

        //(*ddist_dpos)[0]=(*ddist_dpos)[0]*e1.x+(*ddist_dpos)[1]*e2.x+(*ddist_dpos)[2]*e3.x;
        //(*ddist_dpos)[1]=(*ddist_dpos)[0]*e1.y+(*ddist_dpos)[1]*e2.y+(*ddist_dpos)[2]*e3.y;
        //(*ddist_dpos)[2]=(*ddist_dpos)[0]*e1.z+(*ddist_dpos)[1]*e2.z+(*ddist_dpos)[2]*e3.z;

        (*ddist_dp)[offset+0] = tmp.at(0, 0) * X + tmp.at(0, 1) * Y + tmp.at(0, 2) * Z;
        (*ddist_dp)[offset+1] = tmp.at(1, 0) * X + tmp.at(1, 1) * Y + tmp.at(1, 2) * Z;
        (*ddist_dp)[offset+2] = tmp.at(2, 0) * X + tmp.at(2, 1) * Y + tmp.at(2, 2) * Z;
        (*ddist_dp)[offset+3] = tmp.at(3, 0) * X + tmp.at(3, 1) * Y + tmp.at(3, 2) * Z;

        (*ddist_dpos)[0] = X * tmp.at(4, 0) + Y * tmp.at(4, 1) + Z * tmp.at(4, 2);
        (*ddist_dpos)[1] = X * tmp.at(5, 0) + Y * tmp.at(5, 1) + Z * tmp.at(5, 2);
        (*ddist_dpos)[2] = X * tmp.at(6, 0) + Y * tmp.at(6, 1) + Z * tmp.at(6, 2);
      }
      else
      {
        std::vector<float> x;
        x.insert(x.end(), p.begin(), p.end());
        x.push_back(pos.x);
        x.push_back(pos.y);
        x.push_back(pos.z);
        std::vector<float> y(3);
        RotateSdfNode_apply(x.data(), y.data());
        glm::vec3 ps = {y[0], y[1], y[2]};
        d = child->get_distance(ps, ddist_dp, ddist_dpos);
      }

      return d;
    }
    virtual unsigned param_cnt() const override { return 3; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-2*PI,2*PI, ParameterType::DIFFERENTIABLE, "axis_rot_ang_xy"});
      params.push_back({0,-2*PI,2*PI, ParameterType::DIFFERENTIABLE, "axis_rot_ang_z"});
      params.push_back({0,-2*PI,2*PI, ParameterType::DIFFERENTIABLE, "rot_angle"});
      return params;
    }
    virtual AABB get_bbox() const override
    {
      AABB ch_bbox = child->get_bbox();
      std::vector<float> y(3 * 8);
      for (int i = 0; i < 8; ++i)
      {
        std::vector<float> x;
        x.insert(x.end(), p.begin(), p.end());
        if (i % 2) x.push_back(ch_bbox.min_pos.x);
        else x.push_back(ch_bbox.max_pos.x);
        if ((i / 2) % 2) x.push_back(ch_bbox.min_pos.y);
        else x.push_back(ch_bbox.max_pos.y);
        if ((i / 4) % 2) x.push_back(ch_bbox.min_pos.z);
        else x.push_back(ch_bbox.max_pos.z);
        RotateSdfNode_apply(x.data(), y.data()+3*i);
      }
      ch_bbox.min_pos = {y[0], y[1], y[2]};
      ch_bbox.max_pos = {y[0], y[1], y[2]};
      for (int i = 1; i < 8; ++i)
      {
        if (y[i * 3] > ch_bbox.max_pos.x) ch_bbox.max_pos.x = y[i * 3];
        if (y[i * 3] < ch_bbox.min_pos.x) ch_bbox.min_pos.x = y[i * 3];
        if (y[i * 3 + 1] > ch_bbox.max_pos.y) ch_bbox.max_pos.y = y[i * 3 + 1];
        if (y[i * 3 + 1] < ch_bbox.min_pos.y) ch_bbox.min_pos.y = y[i * 3 + 1];
        if (y[i * 3 + 2] > ch_bbox.max_pos.z) ch_bbox.max_pos.z = y[i * 3 + 2];
        if (y[i * 3 + 2] < ch_bbox.min_pos.z) ch_bbox.min_pos.z = y[i * 3 + 2];
      }
      return AABB(ch_bbox.min_pos, ch_bbox.max_pos);
    }
  };

  class OrSdfNode : public TwoChildSdfNode
  {
  public:
    OrSdfNode(unsigned id) : TwoChildSdfNode(id) { name = "Or"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      std::vector<float> ddist_dpos1 = {0,0,0};
      std::vector<float> ddist_dpos2 = {0,0,0};

      int offset_left = ddist_dp ? ddist_dp->size() : 0;
      float d1 = left->get_distance(pos, ddist_dp, &ddist_dpos1);
      int offset_right = ddist_dp ? ddist_dp->size() : 0;
      float d2 = right->get_distance(pos, ddist_dp, &ddist_dpos2);
      int offset_next = ddist_dp ? ddist_dp->size() : 0;

      if (ddist_dp)
      {
        //d(p,x,y) = min(d1(p,x), d2(p,y))
        
        if (d1 < d2)
        {
          for (int i=offset_right;i<offset_next;i++)
            (*ddist_dp)[i] = 0;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = ddist_dpos1[i];
        }
        else
        {
          for (int i=offset_left;i<offset_right;i++)
            (*ddist_dp)[i] = 0;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = ddist_dpos2[i];
        }
      }

      return std::min(d1,d2);
    }
    virtual unsigned param_cnt() const override { return 0; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override { return {}; }
    virtual AABB get_bbox() const override
    {
      AABB b1 = left->get_bbox();
      AABB b2 = right->get_bbox();
      return AABB(glm::min(b1.min_pos, b2.min_pos), glm::max(b1.max_pos, b2.max_pos));
    }
  };

  class AndSdfNode : public TwoChildSdfNode
  {
  public:
    AndSdfNode(unsigned id) : TwoChildSdfNode(id) { name = "And"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      std::vector<float> ddist_dpos1 = {0,0,0};
      std::vector<float> ddist_dpos2 = {0,0,0};

      int offset_left = ddist_dp ? ddist_dp->size() : 0;
      float d1 = left->get_distance(pos, ddist_dp, &ddist_dpos1);
      int offset_right = ddist_dp ? ddist_dp->size() : 0;
      float d2 = right->get_distance(pos, ddist_dp, &ddist_dpos2);
      int offset_next = ddist_dp ? ddist_dp->size() : 0;

      if (ddist_dp)
      {
        //d(p,x,y) = max(d1(p,x), d2(p,y))
        
        if (d1 > d2)
        {
          for (int i=offset_right;i<offset_next;i++)
            (*ddist_dp)[i] = 0;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = ddist_dpos1[i];
        }
        else
        {
          for (int i=offset_left;i<offset_right;i++)
            (*ddist_dp)[i] = 0;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = ddist_dpos2[i];
        }
      }

      return std::max(d1,d2);
    }
    virtual unsigned param_cnt() const override { return 0; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override { return {}; }
    virtual AABB get_bbox() const override
    {
      AABB b1 = left->get_bbox();
      AABB b2 = right->get_bbox();
      return AABB(glm::max(b1.min_pos, b2.min_pos), glm::min(b1.max_pos, b2.max_pos));
    }
  };

  class SubtractSdfNode : public TwoChildSdfNode
  {
  public:
    SubtractSdfNode(unsigned id) : TwoChildSdfNode(id) { name = "Subtract"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      std::vector<float> ddist_dpos1 = {0,0,0};
      std::vector<float> ddist_dpos2 = {0,0,0};

      int offset_left = ddist_dp ? ddist_dp->size() : 0;
      float d1 = left->get_distance(pos, ddist_dp, &ddist_dpos1);
      int offset_right = ddist_dp ? ddist_dp->size() : 0;
      float d2 = right->get_distance(pos, ddist_dp, &ddist_dpos2);
      int offset_next = ddist_dp ? ddist_dp->size() : 0;

      if (ddist_dp)
      {
        //Node1 - Node2
        //d(p,x,y) = max(d1(p,x), -d2(p,y))
        
        if (d1 > -d2)
        {
          for (int i=offset_right;i<offset_next;i++)
            (*ddist_dp)[i] = 0;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = ddist_dpos1[i];
        }
        else
        {
          for (int i=offset_left;i<offset_right;i++)
            (*ddist_dp)[i] = 0;
          for (int i=offset_right;i<offset_next;i++)
            (*ddist_dp)[i] *= -1;
          for (int i=0;i<3;i++)
            (*ddist_dpos)[i] = -ddist_dpos2[i];
        }
      }

      return std::max(d1,-d2);
    }
    virtual unsigned param_cnt() const override { return 0; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override { return {}; }
    virtual AABB get_bbox() const override
    {
      //not always optimal, but valid
      return left->get_bbox();
    }
  };

  ProceduralSdf SdfGenInstance::generate(std::span<const float> parameters)
  {
    for (int i = 0; i < all_params.size(); ++i)
    {
      if (i < parameters.size())
      {
        all_params[i] = parameters[i];
      }
      else
      {
        all_params[i] = 0;
      }
    }
    return ProceduralSdf(root, structure, all_nodes);
  }

  SdfGenInstance::SdfGenInstance(const UPGStructure &structure)
  {
    recreate(structure);
  }
  void SdfGenInstance::recreate(const UPGStructure &_structure)
  {
    structure = _structure;
    all_params.clear();
    all_nodes.clear();
    
    std::vector<SdfNode *> nodes;
    std::vector<std::pair<SdfNode *, unsigned>> param_startings;
    int i = 0;
    do
    {
      unsigned n = 0;
      if (i < structure.s.size())
      {
        n = structure.s[i];
      }
      SdfNode *node = sdf_node_by_node_type_id(n, i);
      all_nodes.push_back(std::unique_ptr<SdfNode>(node));
      desc.add_parameters(node->get_ID(), node->get_node_name(), node->get_parameters_block(scene_bbox));
      param_startings.push_back({node, all_params.size()});
      all_params.resize(all_params.size() + node->param_cnt());
      
      if (i == 0)
      {
        root = node;
        if (node->child_cnt() > 0)
        {
          nodes.push_back(node);
        }
      }
      else
      {
        SdfNode *last = nodes[nodes.size() - 1];
        
        if (!last->add_child(node))
        {
          nodes.pop_back();
        }
        if (node->child_cnt() > 0)
        {
          nodes.push_back(node);
        }
      }
      ++i;
    } while (nodes.size() > 0);

    int offset = 0;
    for (auto &nptr : all_nodes)
    {
      nptr->set_param_span(std::span<my_float>(all_params.data() + offset, nptr->param_cnt()));
      offset += nptr->param_cnt();
    }
  }

  SdfNode *sdf_node_by_node_type_id(uint16_t num, unsigned id)
  {
    SdfNode *node = NULL;
    switch(num)
    {
      case SdfNode::SPHERE: 
        node = new SphereSdfNode(id);
        break;
      case SdfNode::MOVE:
        node = new MoveSdfNode(id);
        break;
      case SdfNode::OR:
        node = new OrSdfNode(id);
        break;
      case SdfNode::BOX:
        node = new BoxSdNode(id);
        break;
      case SdfNode::CYLINDER:
        node = new CylinderSdNode(id);
        break;
      case SdfNode::ROUNDED_BOX:
        node = new RoundBoxSdNode(id);
        break;
      case SdfNode::PRISM:
        node = new Prism(id);
        break;
      case SdfNode::CONE:
        node = new Cone(id);
        break;
      case SdfNode::AND:
        node = new AndSdfNode(id);
        break;
      case SdfNode::SUBTRACT:
        node = new SubtractSdfNode(id);
        break;
      case SdfNode::ROTATE:
        node = new RotateSdfNode(id);
        break;
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;
    }
    return node;
  }

  constexpr int merge_node_num = 3;
  constexpr int move_node_num  = 2;

  int get_node_power(int node_id)
  {
    // how many children each node has. It should be done differently, but now it's just testing
    std::vector<int> node_powers = {0, 0, 1, 2, 0, 0, 0, 0, 0, 2, 2};
    assert(node_id >= 0 && node_id < node_powers.size());
    return node_powers[node_id];
  }

  int get_node_param_count(int node_id)
  {
    // how many parameters each node has. It should be done differently, but now it's just testing
    std::vector<int> node_param_counts = {0, 1, 3, 0, 3, 2, 4, 2, 3, 0, 0};
    assert(node_id >= 0 && node_id < node_param_counts.size());
    return node_param_counts[node_id];
  }

  int get_position_index(const UPGStructure &structure, int start, int p_start)
  {
    int pi = p_start;
    int si = start;
    while (si < structure.s.size() && get_node_power(structure.s[si]) == 1)
    {
      if (structure.s[si] == move_node_num)
        return pi;
      si++;
      pi += get_node_param_count(structure.s[si]);
    }
    return -1;
  }

  int parse_node_rec(const UPGStructure &structure, std::vector<UPGPart> &parts, bool merge_level, int start)
  {

    bool is_merge = merge_level && (structure.s[start] == merge_node_num);
    int power = get_node_power(structure.s[start]);
    int pos = start + 1;

    for (int i = 0; i < power; i++)
    {
      int fin_pos = parse_node_rec(structure, parts, is_merge, pos);
      if (is_merge && structure.s[pos] != merge_node_num)
      {
        int p_st = parts.empty() ? 0 : parts.back().p_range.second;
        int p_cnt = 0;
        for (int p = pos; p < fin_pos; p++)
          p_cnt += get_node_param_count(structure.s[p]);
        parts.push_back(UPGPart({pos, fin_pos}, {p_st, p_st + p_cnt}, get_position_index(structure, pos, p_st)));
      }
      pos = fin_pos;
    }
    return pos;
  }

  int total_param_count(const UPGStructure &structure)
  {
    int cnt = 0;
    for (auto &node_id : structure.s)
      cnt += get_node_param_count(node_id);
    return cnt;
  }

  std::vector<UPGPart> get_sdf_parts(const UPGStructure &structure)
  {
    std::vector<UPGPart> parts;
    parse_node_rec(structure, parts, true, 0);
    if (parts.empty())
      parts.push_back(UPGPart({0, (int)structure.s.size()}, {0, total_param_count(structure)}, get_position_index(structure, 0, 0)));
    for (auto &g : parts)
      logerr("part [%d %d][%d %d] %d", (int)g.s_range.first, (int)g.s_range.second, g.p_range.first, g.p_range.second, g.position_index);
    return parts;
  }
}