#include "sdf_node.h"

namespace upg
{
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

  class SphereSdfNode : public PrimitiveSdfNode
  {
    static constexpr int CENTER_X = 0;
    static constexpr int CENTER_Y = 1;
    static constexpr int CENTER_Z = 2;
    static constexpr int RADIUS = 3;
  public:
    SphereSdfNode(unsigned id) : PrimitiveSdfNode(id) { name = "Sphere"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      float d = std::max(1e-9f, glm::length(glm::vec3(p[CENTER_X], p[CENTER_Y], p[CENTER_Z]) - pos));
      
      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 4);
        (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/d;
        (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/d;
        (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/d;
        (*ddist_dp)[offset+3] = -1;

        (*ddist_dpos)[0] = -(p[CENTER_X] - pos.x)/d;
        (*ddist_dpos)[1] = -(p[CENTER_Y] - pos.y)/d;
        (*ddist_dpos)[2] = -(p[CENTER_Z] - pos.z)/d;
      }

      return d - p[RADIUS];
    }

    virtual unsigned param_cnt() const override { return 4; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_z"});
      params.push_back({1,0.01,10, ParameterType::DIFFERENTIABLE, "radius"});
      return params;
    }
  };

  class BoxSdNode : public PrimitiveSdfNode
  {
    static constexpr int CENTER_X = 0;
    static constexpr int CENTER_Y = 1;
    static constexpr int CENTER_Z = 2;
    
    static constexpr int SIZE_X = 3;
    static constexpr int SIZE_Y = 4;
    static constexpr int SIZE_Z = 5;

  public:
    BoxSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Box"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 size(p[SIZE_X], p[SIZE_Y], p[SIZE_Z]);
      glm::vec3 pos_in_new_space = pos - glm::vec3(p[CENTER_X], p[CENTER_Y], p[CENTER_Z]);
      glm::vec3 q = abs(pos_in_new_space) - size;
      
      float pos_dist = std::max(1e-9f, glm::length(pos_in_new_space));
      float d = glm::length(glm::max(q, glm::vec3(0.0, 0.0, 0.0))) + std::min(std::max(q.y, q.z), 0.f);

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 6);
        (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+3] = -1;
        (*ddist_dp)[offset+4] = -1;
        (*ddist_dp)[offset+5] = -1;

        (*ddist_dpos)[0] = -(p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dpos)[1] = -(p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dpos)[2] = -(p[CENTER_Z] - pos.z)/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 6; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_z"});
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "size_x"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "size_y"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "size_z"});

      return params;
    }
  };

  class CylinderSdNode : public PrimitiveSdfNode
  {
    static constexpr int CENTER_X = 0;
    static constexpr int CENTER_Y = 1;
    static constexpr int CENTER_Z = 2;
    
    static constexpr int HEIGHT = 3;
    static constexpr int RADIUS = 4;

  public:
    CylinderSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Cylinder"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 pos_in_new_space = pos - glm::vec3(p[CENTER_X], p[CENTER_Y], p[CENTER_Z]);
      glm::vec2 vec_d = glm::abs(glm::vec2(glm::length(glm::vec2(pos_in_new_space.x, pos_in_new_space.z)), pos_in_new_space.y)) - glm::vec2(p[RADIUS], p[HEIGHT]);
      
      float pos_dist = std::max(1e-9f, glm::length(pos_in_new_space));
      float d = std::min(std::max(vec_d.x, vec_d.y), 0.f) + glm::length(glm::max(vec_d, glm::vec2(0, 0)));

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 5);
        (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+3] = -1;
        (*ddist_dp)[offset+4] = -1;

        (*ddist_dpos)[0] = -(p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dpos)[1] = -(p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dpos)[2] = -(p[CENTER_Z] - pos.z)/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 5; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_z"});
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});
      params.push_back({2,0.01,10, ParameterType::DIFFERENTIABLE, "radius"});

      return params;
    }
  };

  class RoundedCylinderSdNode : public PrimitiveSdfNode
  {
    static constexpr int CENTER_X = 0;
    static constexpr int CENTER_Y = 1;
    static constexpr int CENTER_Z = 2;
    
    static constexpr int RA = 3;
    static constexpr int RB = 4;
    static constexpr int HEIGHT = 5;

  public:
    RoundedCylinderSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "RoundedCylinder"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 pos_in_new_space = pos - glm::vec3(p[CENTER_X], p[CENTER_Y], p[CENTER_Z]);
      glm::vec2 vec_d = glm::vec2(glm::length(glm::vec2(pos_in_new_space.x, pos_in_new_space.z)) - 2.0 * p[RA] + p[RB], std::abs(pos_in_new_space.y - p[HEIGHT]));
      
      float pos_dist = std::max(1e-9f, glm::length(pos_in_new_space));
      float d = std::min(std::max(vec_d.x, vec_d.y), 0.f) + glm::length(glm::max(vec_d, glm::vec2(0, 0))) - p[RB];

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 6);
        (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+3] = -1;
        (*ddist_dp)[offset+4] = -1;
        (*ddist_dp)[offset+5] = -1;

        (*ddist_dpos)[0] = -(p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dpos)[1] = -(p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dpos)[2] = -(p[CENTER_Z] - pos.z)/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 6; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_z"});
      
      params.push_back({1,0.01,10, ParameterType::DIFFERENTIABLE, "ra"});
      params.push_back({2,0.01,10, ParameterType::DIFFERENTIABLE, "rb"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});

      return params;
    }
  };

  class PyramidSdNode : public PrimitiveSdfNode
  {
    static constexpr int CENTER_X = 0;
    static constexpr int CENTER_Y = 1;
    static constexpr int CENTER_Z = 2;
    
    static constexpr int H = 3;

  public:
    PyramidSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Pyramid"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 pos_in_new_space = pos - glm::vec3(p[CENTER_X], p[CENTER_Y], p[CENTER_Z]);
      float pos_dist = std::max(1e-9f, glm::length(pos_in_new_space));

      float m2 = p[H] * p[H] + 0.25;
      
      pos_in_new_space.x = abs(pos_in_new_space.x);
      pos_in_new_space.z = abs(pos_in_new_space.z);
      pos_in_new_space = (pos_in_new_space.z > pos_in_new_space.x) ? glm::vec3(pos_in_new_space.z, pos_in_new_space.y, pos_in_new_space.x) : pos_in_new_space;
      pos_in_new_space -= glm::vec3(0.5, 0, 0.5);

      glm::vec3 q = glm::vec3(pos_in_new_space.z, 
                              p[H] * pos_in_new_space.y - 0.5 * pos_in_new_space.x, 
                              p[H] * pos_in_new_space.x + 0.5 * pos_in_new_space.y);

      float s = std::max(-q.x, 0.f);
      float t = glm::clamp((q.y - 0.5 * pos_in_new_space.z) / (m2 + 0.25), 0.0, 1.0);

      float a = m2 * (q.x + s) * (q.x + s) + q.y * q.y;
      float b = m2 * (q.x + 0.5 * t) * (q.x + 0.5 * t) + (q.y - m2 * t) * (q.y - m2 * t);
      float d2 = std::min(q.y, -q.x * m2 - q.y * 0.5f) > 0.0 ? 0.0 : std::min(a,b);
      
      float d = sqrt((d2 + q.z * q.z) / m2 ) * glm::sign(std::max(q.z, -pos_in_new_space.y));

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 4);
        (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+3] = -1;

        (*ddist_dpos)[0] = -(p[CENTER_X] - pos.x)/pos_dist;
        (*ddist_dpos)[1] = -(p[CENTER_Y] - pos.y)/pos_dist;
        (*ddist_dpos)[2] = -(p[CENTER_Z] - pos.z)/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 4; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "center_z"});
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});

      return params;
    }
  };

  class MoveSdfNode : public OneChildSdfNode
  {
    static constexpr int MOVE_X = 0;
    static constexpr int MOVE_Y = 1;
    static constexpr int MOVE_Z = 2;
    static constexpr int RADIUS = 3;
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
      float d = child->get_distance(pos - glm::vec3(p[MOVE_X], p[MOVE_Y], p[MOVE_Z]));

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
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_x"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_y"});
      params.push_back({0,-5,5, ParameterType::DIFFERENTIABLE, "move_z"});
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
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const override { return {}; }
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
      desc.add_parameters(node->get_ID(), node->get_node_name(), node->get_parameters_block());
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
        node = new RoundedCylinderSdNode(id);
        break;
      case 7:
        node = new PyramidSdNode(id);
        break;
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;
    }
    return node;
  }
}