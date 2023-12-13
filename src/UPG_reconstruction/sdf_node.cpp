#include "sdf_node.h"
#include <stdio.h>

template<typename T>
T square(T x) { return x * x; }

float __enzyme_autodiffFloat(float (*)(float), float);
double __enzyme_autodiffDouble(double (*)(double), double);
float __enzyme_autodiff(...);

void
f()
{
  printf("float  d/dx %f\n", __enzyme_autodiffFloat(square<float>, 1.0f));
}

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


  float 
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
      float d = std::max(1e-9f, glm::length(pos));
      float p_args[4], d_p[4] = {0};
      p_args[0] = pos.x;
      p_args[1] = pos.y;
      p_args[2] = pos.z;
      p_args[3] = p[RADIUS];

      float ans = __enzyme_autodiff((void*)diff_sphere_sdf, p_args, d_p);

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 1);
        (*ddist_dp)[offset] = d_p[3];
        
        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];
      }

      return d - p[RADIUS];
    }

    virtual unsigned param_cnt() const override { return 1; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      float max_r = 0.5*length(scene_bbox.max_pos-scene_bbox.min_pos);
      std::vector<ParametersDescription::Param> params;
      params.push_back({1,0.01f*max_r,max_r, ParameterType::DIFFERENTIABLE, "radius"});
      return params;
    }
  };

  float
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
      glm::vec3 size(p[SIZE_X], p[SIZE_Y], p[SIZE_Z]);
      glm::vec3 q = abs(pos) - size;
      float d1 = glm::length(glm::max(q, glm::vec3(0.0, 0.0, 0.0)));
      float d2 = std::min(std::max(q.x,std::max(q.y, q.z)), 0.f);
      float d = d1 + d2;

      float p_args[6], d_p[6] = {0};
      p_args[0] = pos.x;
      p_args[1] = pos.y;
      p_args[2] = pos.z;
      p_args[3] = p[SIZE_X];
      p_args[4] = p[SIZE_Y];
      p_args[5] = p[SIZE_Z];

      if (ddist_dp)
      {        
        float ans = __enzyme_autodiff((void*)diff_box_sdf, p_args, d_p);
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 3);

        (*ddist_dp)[offset]   = d_p[3];
        (*ddist_dp)[offset+1] = d_p[4];
        (*ddist_dp)[offset+2] = d_p[5];

        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];

        // std::cout << ans << " " << d << std::endl;
      }

      return d;
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
  };

  float
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
      glm::vec3 size(p[SIZE_X], p[SIZE_Y], p[SIZE_Z]);
      glm::vec3 q = abs(pos) - size;
      float d1 = glm::length(glm::max(q, glm::vec3(0.0, 0.0, 0.0)));
      float d2 = std::min(std::max(q.x,std::max(q.y, q.z)), 0.f);
      float d = d1 + d2 - p[RADIUS];

      float p_args[7], d_p[7] = {0};
      p_args[0] = pos.x;
      p_args[1] = pos.y;
      p_args[2] = pos.z;
      p_args[3] = p[SIZE_X];
      p_args[4] = p[SIZE_Y];
      p_args[5] = p[SIZE_Z];
      p_args[6] = p[RADIUS];

      if (ddist_dp)
      {        
        float ans = __enzyme_autodiff((void*)diff_round_box_sdf, p_args, d_p);
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 4);

        (*ddist_dp)[offset]   = d_p[3];
        (*ddist_dp)[offset+1] = d_p[4];
        (*ddist_dp)[offset+2] = d_p[5];
        (*ddist_dp)[offset+3] = d_p[6];

        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];

        // std::cout << ans << " " << d << std::endl;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 4; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      glm::vec3 size = scene_bbox.max_pos-scene_bbox.min_pos;
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "size_x"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "size_y"});
      params.push_back({5,0.01f*size.z,size.z, ParameterType::DIFFERENTIABLE, "size_z"});
      params.push_back({0.3,0.01f*size.z,size.z * 4, ParameterType::DIFFERENTIABLE, "radius"});

      return params;
    }
  };

  float
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
      glm::vec2 vec_d = glm::abs(glm::vec2(glm::length(glm::vec2(pos.x, pos.z)), pos.y)) - glm::vec2(p[RADIUS], p[HEIGHT]);
      
      float d = std::min(std::max(vec_d.x, vec_d.y), 0.f) + glm::length(glm::max(vec_d, glm::vec2(0, 0)));

      if (ddist_dp)
      {
        float p_args[5], d_p[5] = {0};
        p_args[0] = pos.x;
        p_args[1] = pos.y;
        p_args[2] = pos.z;
        p_args[3] = p[HEIGHT];
        p_args[4] = p[RADIUS];
        
        float ans = __enzyme_autodiff((void*)diff_cylinder_sdf, p_args, d_p);
        ans = diff_cylinder_sdf(p_args);

        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 2);
        
        (*ddist_dp)[offset] = d_p[3];
        (*ddist_dp)[offset+1] = d_p[4];

        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 2; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});
      params.push_back({2,0.01,10, ParameterType::DIFFERENTIABLE, "radius"});

      return params;
    }
  };

  float
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
      glm::vec3 q = abs(pos);
      float d = std::max(q.z - p[H2], std::max(q.x * 0.866025f + pos.y*0.5f, -pos.y) - p[H1] * 0.5f);

      if (ddist_dp)
      {
        float p_args[5], d_p[5] = {0};
        p_args[0] = pos.x;
        p_args[1] = pos.y;
        p_args[2] = pos.z;
        p_args[3] = p[H1];
        p_args[4] = p[H2];

        float ans =  __enzyme_autodiff((void*)diff_prism_sdf, p_args, d_p);

        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 2);
        (*ddist_dp)[offset+0] = d_p[3];
        (*ddist_dp)[offset+1] = d_p[4];

        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 2; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height1"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height2"});

      return params;
    }
  };

  float
  diff_cone_sdf(float params[6])
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

  class Cone : public PrimitiveSdfNode
  {
    static constexpr int SIZE_X = 0;
    static constexpr int SIZE_Y = 1;
    static constexpr int SIZE_Z = 2;

  public:
    Cone(unsigned id) : PrimitiveSdfNode(id) { name = "Cone"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 size(p[SIZE_X], p[SIZE_Y], p[SIZE_Z]);
      glm::vec3 q = abs(pos) - size;
      float d1 = glm::length(glm::max(q, glm::vec3(0.0, 0.0, 0.0)));
      float d2 = std::min(std::max(q.x,std::max(q.y, q.z)), 0.f);
      float d = d1 + d2;

      float p_args[6], d_p[6] = {0};
      p_args[0] = pos.x;
      p_args[1] = pos.y;
      p_args[2] = pos.z;
      p_args[3] = p[SIZE_X];
      p_args[4] = p[SIZE_Y];
      p_args[5] = p[SIZE_Z];

      if (ddist_dp)
      {        
        float ans = __enzyme_autodiff((void*)diff_box_sdf, p_args, d_p);
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 3);

        (*ddist_dp)[offset]   = d_p[3];
        (*ddist_dp)[offset+1] = d_p[4];
        (*ddist_dp)[offset+2] = d_p[5];

        (*ddist_dpos)[0] = d_p[0];
        (*ddist_dpos)[1] = d_p[1];
        (*ddist_dpos)[2] = d_p[2];

        // std::cout << ans << " " << d << std::endl;
      }

      return d;
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
    return ProceduralSdf(root);
  }

  SdfGenInstance::SdfGenInstance(const UPGStructure &structure)
  {
    recreate(structure);
  }
  void SdfGenInstance::recreate(const UPGStructure &structure)
  {
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
      case 1: 
        node = new SphereSdfNode(id);
        break;
      case 2:
        node = new MoveSdfNode(id);
        break;
      case 3:
        node = new OrSdfNode(id);
        break;
      case 4:
        node = new BoxSdNode(id);
        break;
      case 5:
        node = new CylinderSdNode(id);
        break;
      case 6:
        node = new RoundBoxSdNode(id);
        break;
      case 7:
        node = new Prism(id);
        break;
      case 8:
        node = new Cone(id);
        break;
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;
    }
    return node;
  }
}