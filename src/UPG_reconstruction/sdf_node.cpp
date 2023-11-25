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
    all_params.clear();
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
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;

    }
    return node;
  }
}