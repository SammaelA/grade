#include "sdf_node.h"
#include <stdio.h>
#include "generation.h"
#include "autodiff/autodiff.h"
#include "sdf_grid.h"

float __enzyme_autodiff(...);

namespace upg
{
AABB ProceduralSdf::scene_bbox;

#define GET_DISTANCE_WITH_DIFF_BATCH(func, p_cnt)                            \
  {                                                                          \
    float p_args[3+p_cnt];                                                   \
    float d_p[3+p_cnt] = {0};                                                \
    for (int i = 0; i < batch_size; i++)                                     \
    {                                                                        \
      p_args[0] = positions[3 * i + 0];                                      \
      p_args[1] = positions[3 * i + 1];                                      \
      p_args[2] = positions[3 * i + 2];                                      \
      for (int j = 0; j < p_cnt; j++)                                        \
        p_args[3 + j] = p[j];                                                \
      if (ddist_dparams)                                                     \
      {                                                                      \
        __enzyme_autodiff(func, p_args, d_p);                                \
        for (int j = 0; j < p_cnt; j++)                                      \
          ddist_dparams[p_offset * batch_size + p_cnt * i + j] = d_p[3 + j]; \
        ddist_dpos[3 * i + 0] = d_p[0];                                      \
        ddist_dpos[3 * i + 1] = d_p[1];                                      \
        ddist_dpos[3 * i + 2] = d_p[2];                                      \
      }                                                                      \
      distances[i] = func(p_args);                                           \
    }                                                                        \
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      //printf("Ppos %f %f %f\n",positions[3*0+0], positions[3*0+1], positions[3*0+2]);
      for (int i=0;i<batch_size;i++)
        distances[i] = sqrtf(positions[3*i+0]*positions[3*i+0] + positions[3*i+1]*positions[3*i+1] + positions[3*i+2]*positions[3*i+2]) - p[RADIUS];
      if (ddist_dparams)
      {
        for (int i=0;i<batch_size;i++)
          ddist_dparams[p_offset*batch_size + 1*i] = -1;
      }
      if (ddist_dpos)
      {
        for (int i=0;i<batch_size;i++)
        {
          ddist_dpos[3*i+0] = positions[3*i+0]/std::max(1e-9f,distances[i]+p[RADIUS]);
          ddist_dpos[3*i+1] = positions[3*i+1]/std::max(1e-9f,distances[i]+p[RADIUS]);
          ddist_dpos[3*i+2] = positions[3*i+2]/std::max(1e-9f,distances[i]+p[RADIUS]);
        }
      }
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      GET_DISTANCE_WITH_DIFF_BATCH(diff_box_sdf, 3)                                   
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      GET_DISTANCE_WITH_DIFF_BATCH(diff_round_box_sdf, 4)                                   
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      GET_DISTANCE_WITH_DIFF_BATCH(diff_cylinder_sdf, 2)                                   
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      GET_DISTANCE_WITH_DIFF_BATCH(diff_prism_sdf, 2)                                   
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      GET_DISTANCE_WITH_DIFF_BATCH(diff_cone_sdf, 3)                                   
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

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      if (stack.size() - stack_head < 3*batch_size)
        assert(false);
      for (int i=0;i<batch_size;i++)
      {
        stack[stack_head + 3*i+0] = positions[3*i+0] - p[MOVE_X];
        stack[stack_head + 3*i+1] = positions[3*i+1] - p[MOVE_Y];
        stack[stack_head + 3*i+2] = positions[3*i+2] - p[MOVE_Z];
      }
      child->get_distance_batch(batch_size, stack.data() + stack_head, distances, ddist_dparams, ddist_dpos, stack, stack_head + 3*batch_size);
      if (ddist_dparams)
      {
        for (int i=0;i<batch_size;i++)
        {
          ddist_dparams[p_offset*batch_size + 3*i+0] = -ddist_dpos[3*i+0];
          ddist_dparams[p_offset*batch_size + 3*i+1] = -ddist_dpos[3*i+1];
          ddist_dparams[p_offset*batch_size + 3*i+2] = -ddist_dpos[3*i+2];
        }
      }
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

  void RotateSdfNode_apply(const float *in, float *out)
  {
    float x = cosf(in[0]) * cosf(in[1]);
    float y = sinf(in[0]) * cosf(in[1]);
    float z = sinf(in[1]);
    float angle = in[2];
    float pos_x = in[3];
    float pos_y = in[4];
    float pos_z = in[5];
    float c = cosf(angle);
    float s = sinf(angle);

    out[0] = (c + (1 - c) * x * x)     * pos_x + ((1 - c) * x * y + s * z) * pos_y + ((1 - c) * x * z - s * y) * pos_z;
    out[1] = ((1 - c) * x * y - s * z) * pos_x + (c + (1 - c) * y * y)     * pos_y + ((1 - c) * y * z + s * x) * pos_z;
    out[2] = ((1 - c) * x * z + s * y) * pos_x + ((1 - c) * z * y - s * x) * pos_y + (c + (1 - c) * z * z)     * pos_z;
  }

  void get_rotate_mat(const float *in, float *out)
  {
    float x = cosf(in[0]) * cosf(in[1]);
    float y = sinf(in[0]) * cosf(in[1]);
    float z = sinf(in[1]);
    float c = cosf(in[2]);
    float s = sinf(in[2]);
    out[0] = (c + (1 - c) * x * x);
    out[1] = ((1 - c) * x * y + s * z);
    out[2] = ((1 - c) * x * z - s * y);
    out[3] = ((1 - c) * x * y - s * z);
    out[4] = (c + (1 - c) * y * y);
    out[5] = ((1 - c) * y * z + s * x);
    out[6] = ((1 - c) * x * z + s * y);
    out[7] = ((1 - c) * z * y - s * x);
    out[8] = (c + (1 - c) * z * z);
  }

  class RotateSdfNode : public OneChildSdfNode
  {
    static constexpr int AXIS_ANG_XY = 0;
    static constexpr int AXIS_ANG_Z = 1;
    static constexpr int ANGLE = 2;
  public:
    RotateSdfNode(unsigned id) : OneChildSdfNode(id) { name = "Rotate"; }

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      float rot[9];
      float rot_jac[27];
      if (ddist_dparams)
      {
        ENZYME_EVALUATE_WITH_DIFF(get_rotate_mat, 3, 9, p.data(), rot, rot_jac);
      }
      else
        get_rotate_mat(p.data(), rot);
      
      if (stack.size() - stack_head < 3*batch_size)
        assert(false);
      for (int i=0;i<batch_size;i++)
      {
        stack[stack_head + 3*i+0] = rot[3*0+0]*positions[3*i+0] + rot[3*0+1]*positions[3*i+1] + rot[3*0+2]*positions[3*i+2];
        stack[stack_head + 3*i+1] = rot[3*1+0]*positions[3*i+0] + rot[3*1+1]*positions[3*i+1] + rot[3*1+2]*positions[3*i+2];
        stack[stack_head + 3*i+2] = rot[3*2+0]*positions[3*i+0] + rot[3*2+1]*positions[3*i+1] + rot[3*2+2]*positions[3*i+2];
      }

      child->get_distance_batch(batch_size, stack.data() + stack_head, distances, ddist_dparams, ddist_dpos, stack, stack_head + 3*batch_size);
      
      if (ddist_dparams)
      {
        for (int i=0;i<batch_size;i++)
        {
          float darg_dp[3*3];
          for (int j=0;j<9;j++)
            darg_dp[j] = rot_jac[3*j+0]*positions[3*i+0] + rot_jac[3*j+1]*positions[3*i+1] + rot_jac[3*j+2]*positions[3*i+2];
          
          for (int j=0;j<3;j++)
            ddist_dparams[p_offset*batch_size + 3*i+j] = ddist_dpos[3*i+0]*darg_dp[3*j+0] + ddist_dpos[3*i+1]*darg_dp[3*j+1] + ddist_dpos[3*i+2]*darg_dp[3*j+2];

          float dp[3] = {ddist_dpos[3*i+0], ddist_dpos[3*i+1], ddist_dpos[3*i+2]};

          for (int j=0;j<3;j++)
            ddist_dpos[3*i+j] = dp[0]*rot[3*0+j] + dp[1]*rot[3*1+j] + dp[2]*rot[3*2+j];                             
        }
      }
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
        if (i % 2 == 0) x.push_back(ch_bbox.min_pos.x);
        else x.push_back(ch_bbox.max_pos.x);
        if ((i / 2) % 2 == 0) x.push_back(ch_bbox.min_pos.y);
        else x.push_back(ch_bbox.max_pos.y);
        if ((i / 4) % 2 == 0) x.push_back(ch_bbox.min_pos.z);
        else x.push_back(ch_bbox.max_pos.z);
        RotateSdfNode_apply(x.data(), &y.data()[3*i]);
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
    
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      unsigned d1_head = stack_head;
      unsigned d2_head = d1_head + batch_size;
      unsigned pos1_head = d2_head + batch_size;
      unsigned pos2_head = pos1_head + (ddist_dpos!=nullptr)*3*batch_size;
      unsigned child_head = pos2_head + (ddist_dpos!=nullptr)*3*batch_size;

      if (stack.size() < child_head)
        assert(false);

      left->get_distance_batch(batch_size, positions, stack.data() + d1_head, ddist_dparams, 
                               stack.data() + pos1_head, stack, child_head);                            

      right->get_distance_batch(batch_size, positions, stack.data() + d2_head, ddist_dparams, 
                                stack.data() + pos2_head, stack, child_head);
    
      for (int i=0;i<batch_size;i++)
      {
        distances[i] = std::min(stack[d1_head + i], stack[d2_head + i]);
      }
      
      if (ddist_dparams)
      {
        unsigned pcnt_1 = left->subgraph_param_cnt;
        unsigned pcnt_2 = right->subgraph_param_cnt;

        for (int i=0;i<batch_size;i++)
        {
          if (stack[d1_head + i] < stack[d2_head + i])
          {
            for (auto *cn : right->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
          }
          else
          {
            for (auto *cn : left->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
          }
        }
      }
      if (ddist_dpos)
      {
        for (int i=0;i<batch_size;i++)
          for (int j=0;j<3;j++)
            ddist_dpos[3*i + j] = stack[d1_head + i] < stack[d2_head + i] ? stack[pos1_head + 3*i+j] : stack[pos2_head + 3*i+j];
      }

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

    
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      unsigned d1_head = stack_head;
      unsigned d2_head = d1_head + batch_size;
      unsigned pos1_head = d2_head + batch_size;
      unsigned pos2_head = pos1_head + (ddist_dpos!=nullptr)*3*batch_size;
      unsigned child_head = pos2_head + (ddist_dpos!=nullptr)*3*batch_size;

      if (stack.size() < child_head)
        assert(false);

      left->get_distance_batch(batch_size, positions, stack.data() + d1_head, ddist_dparams, 
                               stack.data() + pos1_head, stack, child_head);                            

      right->get_distance_batch(batch_size, positions, stack.data() + d2_head, ddist_dparams, 
                                stack.data() + pos2_head, stack, child_head);
    
      for (int i=0;i<batch_size;i++)
      {
        distances[i] = std::max(stack[d1_head + i], stack[d2_head + i]);
      }
      
      if (ddist_dparams)
      {
        unsigned pcnt_1 = left->subgraph_param_cnt;
        unsigned pcnt_2 = right->subgraph_param_cnt;

        for (int i=0;i<batch_size;i++)
        {
          if (stack[d1_head + i] > stack[d2_head + i])
          {
            for (auto *cn : right->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
          }
          else
          {
            for (auto *cn : left->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
          }
        }
      }
      if (ddist_dpos)
      {
        for (int i=0;i<batch_size;i++)
          for (int j=0;j<3;j++)
            ddist_dpos[3*i + j] = stack[d1_head + i] > stack[d2_head + i] ? stack[pos1_head + 3*i+j] : stack[pos2_head + 3*i+j];
      }

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
    
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      unsigned d1_head = stack_head;
      unsigned d2_head = d1_head + batch_size;
      unsigned pos1_head = d2_head + batch_size;
      unsigned pos2_head = pos1_head + (ddist_dpos!=nullptr)*3*batch_size;
      unsigned child_head = pos2_head + (ddist_dpos!=nullptr)*3*batch_size;

      if (stack.size() < child_head)
        assert(false);

      left->get_distance_batch(batch_size, positions, stack.data() + d1_head, ddist_dparams, 
                               stack.data() + pos1_head, stack, child_head);                            

      right->get_distance_batch(batch_size, positions, stack.data() + d2_head, ddist_dparams, 
                                stack.data() + pos2_head, stack, child_head);
    
      for (int i=0;i<batch_size;i++)
      {
        distances[i] = std::max(stack[d1_head + i], -stack[d2_head + i]);
      }
      
      if (ddist_dparams)
      {
        unsigned pcnt_1 = left->subgraph_param_cnt;
        unsigned pcnt_2 = right->subgraph_param_cnt;

        for (int i=0;i<batch_size;i++)
        {
          if (stack[d1_head + i] > -stack[d2_head + i])
          {
            for (auto *cn : right->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
          }
          else
          {
            for (auto *cn : left->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] = 0;
            }
            for (auto *cn : right->subgraph)
            {
              for (int j=0;j<cn->param_cnt();j++)
                ddist_dparams[batch_size*cn->p_offset + i*cn->param_cnt() + j] *= -1;
            }
          }
        }
      }
      if (ddist_dpos)
      {
        for (int i=0;i<batch_size;i++)
          for (int j=0;j<3;j++)
            ddist_dpos[3*i + j] = stack[d1_head + i] > -stack[d2_head + i] ? stack[pos1_head + 3*i+j] : -stack[pos2_head + 3*i+j];
      }

    }

    virtual unsigned param_cnt() const override { return 0; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override { return {}; }
    virtual AABB get_bbox() const override
    {
      //not always optimal, but valid
      return left->get_bbox();
    }
  };

  class ScaleSdfNode : public OneChildSdfNode
  {
    static constexpr int SCALE = 0;
  public:
    ScaleSdfNode(unsigned id) : OneChildSdfNode(id) { name = "Scale"; }

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      if (stack.size() - stack_head < 3*batch_size)
        assert(false);
      float si = 1/p[SCALE];
      for (int i=0;i<batch_size;i++)
      {
        stack[stack_head + 3*i+0] = si*positions[3*i+0];
        stack[stack_head + 3*i+1] = si*positions[3*i+1];
        stack[stack_head + 3*i+2] = si*positions[3*i+2];
      }
      child->get_distance_batch(batch_size, stack.data() + stack_head, distances, ddist_dparams, ddist_dpos, stack, stack_head + 3*batch_size);
      
      if (ddist_dparams)
      {
        for (int i=0;i<batch_size;i++)
        {
          ddist_dparams[p_offset*batch_size + 3*i+0] = distances[i] - si*(positions[3*i+0]*ddist_dpos[3*i+0] + 
            positions[3*i+1]*ddist_dpos[3*i+1] + positions[3*i+2]*ddist_dpos[3*i+2]);
        }
      }

      for (int i=0;i<batch_size;i++)
        distances[i] *= p[SCALE];
    }

    virtual unsigned param_cnt() const override { return 1; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({1,0.01, 10, ParameterType::DIFFERENTIABLE, "scale"});
      return params;
    }
    virtual AABB get_bbox() const override
    {
      AABB ch_bbox = child->get_bbox();
      return AABB(p[SCALE]*ch_bbox.min_pos, p[SCALE]*ch_bbox.max_pos);
    }
  };

  void ProceduralSdf::set_parameters(std::span<const float> parameters)
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
  }

  void ProceduralSdf::get_distance_batch(unsigned batch_size, float *const positions, float *distances,
                                         float * ddist_dparams, float * ddist_dpos) const
  {
    if (ddist_dparams_transp.size() < batch_size*all_params.size())
      ddist_dparams_transp.resize(batch_size*all_params.size());
    
    unsigned max_stack_size = batch_size*256*(1 + log2(all_nodes.size())); //TODO: calculate it somehow else
    if (stack.size() < max_stack_size)
      stack.resize(max_stack_size, 0.0f);

    root->get_distance_batch(batch_size, positions, distances, ddist_dparams ? ddist_dparams_transp.data() : nullptr, ddist_dpos, stack, 0);
    if (ddist_dparams)
    {
      for (auto &n : all_nodes)
      {
        unsigned p_cnt = n->param_cnt();
        unsigned p_off = n->p_offset;
        if (p_cnt > 0)
        {
          for (int i=0;i<batch_size;i++)
            for (int j=0;j<p_cnt;j++)
            {
              //printf("%d %d %u %d %u\n",i,j,batch_size,(int)all_params.size(), p_off);
              ddist_dparams[i*all_params.size() + p_off + j] = ddist_dparams_transp[batch_size*p_off + i*p_cnt + j];
            }
        }
      }
    }
  }

  float ProceduralSdf::get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp, 
                                    std::vector<float> *ddist_dpos) const
  {
    //return root->get_distance(pos, ddist_dp, ddist_dpos);
    unsigned batch_size = 1;
    if (ddist_dparams_transp.size() < batch_size*all_params.size())
      ddist_dparams_transp.resize(batch_size*all_params.size());
    
    unsigned max_stack_size = batch_size*64*(1 + log2(all_nodes.size())); //TODO: calculate it somehow else
    if (stack.size() < max_stack_size)
      stack.resize(max_stack_size, 0.0f);

    float d = 1e7;
    root->get_distance_batch(1, (float*)(&pos), &d, ddist_dp ?   ddist_dparams_transp.data() : nullptr, 
                                                    ddist_dpos ? ddist_dpos->data() : nullptr, stack, 0);

    if (ddist_dp)
    {
      ddist_dp->resize(batch_size*all_params.size());
      for (auto &n : all_nodes)
      {
        unsigned p_cnt = n->param_cnt();
        unsigned p_off = n->p_offset;
        if (p_cnt > 0)
        {
          for (int i=0;i<batch_size;i++)
            for (int j=0;j<p_cnt;j++)
              (*ddist_dp)[i*all_params.size() + p_off + j] = ddist_dparams_transp[batch_size*p_off + i*p_cnt + j];
        }
      }
    }
    //float dorg = root->get_distance(pos, nullptr, nullptr);
    //if (d != dorg)
    //printf("d dorg %f %f %lu (%f %f %f)\n",d,dorg,(uint64_t)ddist_dp,pos.x, pos.y, pos.z);
    return d;
  }

  ProceduralSdf::ProceduralSdf(const UPGStructure &structure)
  {
    recreate(structure);
  }

  ProceduralSdf::ProceduralSdf(const ProceduralSdf &sdf)
  {
    recreate(sdf.structure);
    set_parameters(sdf.all_params);
  }

  void set_subgraph_params_cnt_rec(SdfNode *node)
  {
    unsigned cnt = node->param_cnt();
    for (auto *c : node->get_children())
    {
      set_subgraph_params_cnt_rec((SdfNode *)c);
      node->subgraph.insert(node->subgraph.end(), c->subgraph.begin(), c->subgraph.end());
      cnt += c->subgraph_param_cnt;
    }
    node->subgraph_param_cnt = cnt;
    node->subgraph.push_back(node);
  }

  void ProceduralSdf::recreate(const UPGStructure &_structure)
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
      nptr->set_param_span(std::span<float>(all_params.data() + offset, nptr->param_cnt()), offset);
      offset += nptr->param_cnt();
    }

    set_subgraph_params_cnt_rec(root);
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
      case SdfNode::GRID:
        node = new GridSdfNode(id, 32, AABB({-1,-1,-1},{1,1,1}));
        break;
      case SdfNode::SCALE:
        node = new ScaleSdfNode(id);
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
    std::vector<int> node_powers = {0, 0, 1, 2, 0, 0, 0, 0, 0, 2, 2, 1};
    assert(node_id >= 0 && node_id < node_powers.size());
    return node_powers[node_id];
  }

  int get_node_param_count(int node_id)
  {
    // how many parameters each node has. It should be done differently, but now it's just testing
    std::vector<int> node_param_counts = {0, 1, 3, 0, 3, 2, 4, 2, 3, 0, 0, 3};
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
    //for (auto &g : parts)
    //  logerr("part [%d %d][%d %d] %d", (int)g.s_range.first, (int)g.s_range.second, g.p_range.first, g.p_range.second, g.position_index);
    return parts;
  }
}