#include "sdf_node.h"

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

  class SphereSdfNode : public PrimitiveSdfNode
  {
    static constexpr int RADIUS = 0;
  public:
    SphereSdfNode(unsigned id) : PrimitiveSdfNode(id) { name = "Sphere"; }
    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      float d = std::max(1e-9f, glm::length(pos));
      
      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 1);
        (*ddist_dp)[offset] = -1;

        (*ddist_dpos)[0] = pos.x/d;
        (*ddist_dpos)[1] = pos.y/d;
        (*ddist_dpos)[2] = pos.z/d;
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

      if (ddist_dp)
      {        
        float dd1_dqx = q.x < 0 ? 0 : q.x/std::max(1e-9f,d1);
        float dd1_dqy = q.y < 0 ? 0 : q.y/std::max(1e-9f,d1);
        float dd1_dqz = q.z < 0 ? 0 : q.z/std::max(1e-9f,d1);

        float dd2_dqx = (q.x<=0 && q.x>=q.y && q.x>=q.z);
        float dd2_dqy = (q.y<=0 && q.y>=q.x && q.y>=q.z);
        float dd2_dqz = (q.z<=0 && q.z>=q.x && q.z>=q.y);

        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 3);
        (*ddist_dp)[offset]   = -(dd1_dqx+dd2_dqx);
        (*ddist_dp)[offset+1] = -(dd1_dqy+dd2_dqy);
        (*ddist_dp)[offset+2] = -(dd1_dqz+dd2_dqz);

        (*ddist_dpos)[0] = glm::sign(pos.x)*(dd1_dqx+dd2_dqx);
        (*ddist_dpos)[1] = glm::sign(pos.y)*(dd1_dqy+dd2_dqy);
        (*ddist_dpos)[2] = glm::sign(pos.z)*(dd1_dqz+dd2_dqz);
                
        //debug("%f %f %f %f - %f %f %f - %f %f %f\n", d, pos.x, pos.y, pos.z, (*ddist_dp)[offset], (*ddist_dp)[offset+1],
        //        (*ddist_dp)[offset+2], (*ddist_dpos)[0], (*ddist_dpos)[1], (*ddist_dpos)[2]);

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
      
      float pos_dist = std::max(1e-9f, glm::length(pos));
      float d = std::min(std::max(vec_d.x, vec_d.y), 0.f) + glm::length(glm::max(vec_d, glm::vec2(0, 0)));

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 2);
        
        float d1 = glm::length(glm::vec2(pos.x, pos.z)) - p[RADIUS];
        float d2 = abs(pos.y) - p[HEIGHT];
        float l2 = glm::length(glm::max(vec_d, glm::vec2(0, 0)));

        //? Calculate derivatives of cylinder radius and height parameters 
        float L1_h = (d1 <= d2 && d2 < 0) ? -1 : 0;
        float L2_h = (d2 > 0) ? -d2 / l2 : 0;

        (*ddist_dp)[offset] = L1_h + L2_h;

        float L1_r = (d1 > d2 && d1 < 0) ? -1 : 0;
        float L2_r = (d1 > 0) ? -d1 / l2 : 0;

        (*ddist_dp)[offset+1] = L1_r + L2_r;

        //? Calculate derivatives of point position parameters
        float const_val = pow(std::max(d2, 0.f), 2);
        float L1_x = (d1 > d2 && d1 < 0) ? pos.x / (d1 + p[RADIUS]) : 0;
        float L2_x = (d1 > 0) ? pos.x * d1 / ((d1 + p[RADIUS]) * sqrt(pow(d1 + p[RADIUS], 2) - 2 * p[RADIUS] * (d1 + p[RADIUS]) + pow(p[RADIUS], 2) + const_val)) : 0;
        
        (*ddist_dpos)[0] = L1_x + L2_x;

        const_val = pow(std::max(d1, 0.f), 2);
        float L1_y = (d1 <= d2 && d2 < 0) ? pos.y / abs(pos.y) : 0;
        float L2_y = (d2 > 0) ? pos.y / abs(pos.y) * d2 / sqrt(const_val + pow(d2, 2)) : 0;
        
        (*ddist_dpos)[1] = L1_y + L2_y;

        const_val = pow(std::max(d2, 0.f), 2);
        float L1_z = (d1 > d2 && d1 < 0) ? pos.z / (d1 + p[RADIUS]) : 0;
        float L2_z = (d1 > 0) ? pos.z * d1 / ((d1 + p[RADIUS]) * sqrt(pow(d1 + p[RADIUS], 2) - 2 * p[RADIUS] * (d1 + p[RADIUS]) + pow(p[RADIUS], 2) + const_val)) : 0;
        
        (*ddist_dpos)[2] = L1_z + L2_z;
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

  class RoundedCylinderSdNode : public PrimitiveSdfNode
  {    
    static constexpr int RA = 0;
    static constexpr int RB = 1;
    static constexpr int HEIGHT = 2;

  public:
    RoundedCylinderSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "RoundedCylinder"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec2 vec_d = glm::vec2(glm::length(glm::vec2(pos.x, pos.z)) - 2.0 * p[RA] + p[RB], std::abs(pos.y - p[HEIGHT]));
      
      float pos_dist = std::max(1e-9f, glm::length(pos));
      float d = std::min(std::max(vec_d.x, vec_d.y), 0.f) + glm::length(glm::max(vec_d, glm::vec2(0, 0))) - p[RB];

      if (ddist_dp)
      {
        int offset = ddist_dp->size();
        ddist_dp->resize(offset + 3);
        // (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        // (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        // (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+0] = -1;
        (*ddist_dp)[offset+1] = -1;
        (*ddist_dp)[offset+2] = -1;

        (*ddist_dpos)[0] = pos.x/pos_dist;
        (*ddist_dpos)[1] = pos.y/pos_dist;
        (*ddist_dpos)[2] = pos.z/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 3; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({1,0.01,10, ParameterType::DIFFERENTIABLE, "ra"});
      params.push_back({2,0.01,10, ParameterType::DIFFERENTIABLE, "rb"});
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});

      return params;
    }
  };

  class PyramidSdNode : public PrimitiveSdfNode
  {
    static constexpr int H = 0;

  public:
    PyramidSdNode(unsigned id) : PrimitiveSdfNode(id) { name = "Pyramid"; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const override
    {
      glm::vec3 pos_in_new_space = pos;
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
        ddist_dp->resize(offset + 1);
        // (*ddist_dp)[offset+0] = (p[CENTER_X] - pos.x)/pos_dist;
        // (*ddist_dp)[offset+1] = (p[CENTER_Y] - pos.y)/pos_dist;
        // (*ddist_dp)[offset+2] = (p[CENTER_Z] - pos.z)/pos_dist;
        (*ddist_dp)[offset+0] = -1;

        (*ddist_dpos)[0] = pos.x/pos_dist;
        (*ddist_dpos)[1] = pos.y/pos_dist;
        (*ddist_dpos)[2] = pos.z/pos_dist;
      }

      return d;
    }

    virtual unsigned param_cnt() const override { return 1; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      params.push_back({5,0.01,10, ParameterType::DIFFERENTIABLE, "height"});

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
        node = new RoundedCylinderSdNode(id);
        break;
      case 7:
        node = new PyramidSdNode(id);
        break;
      case 8:
        node = new AndSdfNode(id);
        break;
      case 9:
        node = new SubtractSdfNode(id);
        break;
      default:
        logerr("invalid node_id %u\n",id);
        node = nullptr;
        break;
    }
    return node;
  }
}