#include "sdf_node.h"
#include <stdio.h>
#include "generation.h"
#include "param_node.h"
#include "autodiff/autodiff.h"
#include "sdf_grid.h"

float __enzyme_autodiff(...);

namespace upg
{
  AABB ProceduralSdf::scene_bbox = AABB({-1,-1,-1},{1,1,1});
  extern std::vector<SdfNodeProperties> node_properties;


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
      for (int j = 0; j < 3+p_cnt; j++)                                      \
        d_p[j] = 0;                                                          \
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

  class OneChildSdfNode : public SdfNode
  {
  protected:
    SdfNode *child;
  public:
    OneChildSdfNode() : SdfNode() { child = NULL; }
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
    TwoChildSdfNode() : SdfNode() { left = NULL; right = NULL; }
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

  class NullSdfNode : public OneChildSdfNode
  {
  protected:
    mutable int childs_params_start = 0, childs_params_end = 0;
    float *global_ddist_dparams = nullptr;
  public:
    NullSdfNode() : OneChildSdfNode() { name = "Null"; }

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      childs_params_start = p_offset;
      child->get_distance_batch(batch_size, positions, distances, global_ddist_dparams, ddist_dpos, stack, stack_head);
      childs_params_end = p_offset + subgraph_param_cnt;
    }

    std::pair <int, int> get_child_param_idxs()
    {
      return std::pair<int, int> {childs_params_start, childs_params_end};
    } 
    void set_global_ddist_dparams(float *ptr)
    {
      global_ddist_dparams = ptr;
    }
    virtual unsigned param_cnt() const override { return 0; }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      return {};
    }
    virtual AABB get_bbox() const override
    {
      return child->get_bbox();
    }
  };

  class AbstractComplexSdfNode : public SdfNode
  {
  protected:
    std::vector<const SdfNode *> childs;
    std::vector<NullSdfNode *> null_nodes;
    std::vector<std::unique_ptr<SdfNode>> inside_nodes;
    SdfNode *root;
    std::vector<float> all_params;
    ParamGenInstance param_inst;
    unsigned input_size;
    unsigned output_param_size() const
    {
      unsigned res = 0;
      for (int i = 0; i < structure.s.size(); ++i)
      {
        unsigned n = 0;
        if (i < structure.s.size())
        {
          n = structure.s[i];
        }
        if (n != 0)
        {
          SdfNode *node = create_node(n);
          res += node->param_cnt();
          delete node;
        }
      }
      return res;
    }
    UPGStructure structure;
    std::vector<std::vector<float>> param_structure;
    SdfNode *create_node_inside_complex_node(uint16_t num)
    {
      if (num == SdfNodeType::UNDEFINED)
      {
        NullSdfNode *node = NULL;
        node = new NullSdfNode();
        null_nodes.push_back(node);
        return node;
      }
      return create_node(num);
    }

    void create(std::vector<std::vector<float>> p_s, UPGStructure s, unsigned in_size)
    {
      param_inst.recreate(p_s, in_size, output_param_size());
      all_params.clear();
      inside_nodes.clear();
      std::vector<SdfNode *> nodes;
      std::vector<std::pair<SdfNode *, unsigned>> param_startings;
      int i = 0;
      do
      {
        unsigned n = 0;
        if (i < s.s.size())
        {
          n = s.s[i];
        }
        SdfNode *node = create_node_inside_complex_node(n);
        inside_nodes.push_back(std::unique_ptr<SdfNode>(node));
        param_startings.push_back({node, all_params.size()});
        all_params.resize(all_params.size() + node->param_cnt());
        if (i == 0)
        {
          root = node;
          if (node->child_cnt() > 0 && n != SdfNodeType::UNDEFINED)
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
          if (node->child_cnt() > 0 && n != SdfNodeType::UNDEFINED)
          {
            nodes.push_back(node);
          }
        }
        ++i;
      } while (nodes.size() > 0);
      int offset = 0;
      for (auto &nptr : inside_nodes)
      {
        nptr->set_param_span(std::span<float>(all_params.data() + offset, nptr->param_cnt()), offset);
        offset += nptr->param_cnt();
      }
    }


  public:
    AbstractComplexSdfNode(std::vector<std::vector<float>> p_s, UPGStructure s, unsigned in_size) : 
      SdfNode(), 
      param_inst({}, 0, 0)
    {
      childs = {};
      param_structure = p_s;
      structure = s;
      input_size = in_size;
      create(p_s, s, in_size);
    }
    bool add_child(SdfNode *node) override 
    {
      if (child_cnt() > childs.size())
        {
          null_nodes[childs.size()]->add_child(node);
          childs.push_back(node);
        }
        return child_cnt() > childs.size();
    }
    std::vector<const SdfNode *> get_children() const override
    {
      return childs;
    }
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      for (auto *n : null_nodes)
        n->set_global_ddist_dparams(ddist_dparams);
/*
      ParamGenInstance inst(param_structure, input_size, output_param_size());
      ParamsGraph param_graph = inst.get_graph(p);
      std::vector<float> pars, jac;
      if (ddist_dpos != nullptr) pars = param_graph.get_params(&jac);
      else pars = param_graph.get_params();

      std::vector<int> param_pairs_idxs = {0};
      if (ddist_dpos != nullptr)
      {
        param_pairs_idxs.push_back(ddist_dpos->size());
      }
      float distance = root->get_distance(pos, ddist_dp, ddist_dpos);
      if (ddist_dpos != nullptr)
      {
        for (int i = 0; i < null_nodes.size(); ++i)
        {
          auto tmp = null_nodes[i]->get_child_param_idxs();
          param_pairs_idxs.push_back(tmp.first);
          param_pairs_idxs.push_back(tmp.second);
        }
      }
      if (ddist_dpos != nullptr)
      {
        param_pairs_idxs.push_back(ddist_dpos->size());
      }
      std::vector<float> buf(0), inside(0);
      if (ddist_dpos != nullptr)
      {
        for (int i = 0; i < param_pairs_idxs.size() - 1; ++i)
        {
          if (i % 2)
          {
            for (int j = param_pairs_idxs[i]; j < param_pairs_idxs[i + 1]; ++j)
            {
              buf.push_back((*ddist_dpos)[j]);
            }
          }
          else
          {
            for (int j = param_pairs_idxs[i]; j < param_pairs_idxs[i + 1]; ++j)
            {
              inside.push_back((*ddist_dpos)[j]);
            }
          }
        }
        
        std::vector <float> new_ddist_part(jac.size() / inside.size());
        int new_size = jac.size() / inside.size();
        for (int i = 0; i < inside.size(); ++i)
        {
          for (int j = 0; j < new_size; ++j)
          {
            if (i == 0)
            {
              new_ddist_part[j] = 0;
            }
            new_ddist_part[j] += jac[i * new_size + j] * inside[i];
          }
        }
        buf.insert(buf.end(), new_ddist_part.begin(), new_ddist_part.end());
        ddist_dpos->clear();
        for (int i = 0; i < buf.size(); ++i)
        {
          ddist_dpos->push_back(buf[i]);
      }
      return distance;
*/
    }
    virtual AABB get_bbox() const override
    {
      return root->get_bbox();
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
    SphereSdfNode() : PrimitiveSdfNode() { name = "Sphere"; }

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

  class BoxSdfNode : public PrimitiveSdfNode
  {
    static constexpr int SIZE_X = 0;
    static constexpr int SIZE_Y = 1;
    static constexpr int SIZE_Z = 2;

  public:
    BoxSdfNode() : PrimitiveSdfNode() { name = "Box"; }

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

  class RoundBoxSdfNode : public PrimitiveSdfNode
  {
    static constexpr int SIZE_X = 0;
    static constexpr int SIZE_Y = 1;
    static constexpr int SIZE_Z = 2;
    static constexpr int RADIUS = 3;

  public:
    RoundBoxSdfNode() : PrimitiveSdfNode() { name = "RoundBox"; }

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

  class CylinderSdfNode : public PrimitiveSdfNode
  {
    static constexpr int HEIGHT = 0;
    static constexpr int RADIUS = 1;

  public:
    CylinderSdfNode() : PrimitiveSdfNode() { name = "Cylinder"; }

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

  class PrismSdfNode : public PrimitiveSdfNode
  {
    static constexpr int H1 = 0;
    static constexpr int H2 = 1;

  public:
    PrismSdfNode() : PrimitiveSdfNode() { name = "Prism"; }

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

  class ConeSdfNode : public PrimitiveSdfNode
  {
    static constexpr int C1 = 0;
    static constexpr int C2 = 1;
    static constexpr int HEIGHT = 2;

  public:
    ConeSdfNode() : PrimitiveSdfNode() { name = "Cone"; }

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
    MoveSdfNode() : OneChildSdfNode() { name = "Move"; }

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
        unsigned cnt = param_cnt();
        for (int i=0;i<batch_size;i++)
        {
          ddist_dparams[p_offset*batch_size + cnt*i+0] = -ddist_dpos[cnt*i+0];
          ddist_dparams[p_offset*batch_size + cnt*i+1] = -ddist_dpos[cnt*i+1];
          ddist_dparams[p_offset*batch_size + cnt*i+2] = -ddist_dpos[cnt*i+2];
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
    RotateSdfNode() : OneChildSdfNode() { name = "Rotate"; }

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
    OrSdfNode() : TwoChildSdfNode() { name = "Or"; }
    
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
    AndSdfNode() : TwoChildSdfNode() { name = "And"; }

    
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
    SubtractSdfNode() : TwoChildSdfNode() { name = "Subtract"; }
    
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
    ScaleSdfNode() : OneChildSdfNode() { name = "Scale"; }

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

  //test start
  class ChairSdfNode : public AbstractComplexSdfNode
  {
  public:
    ChairSdfNode() : AbstractComplexSdfNode(param_structure(), structure(), 6) { name = "chair"; }//last param is param_cnt
    virtual unsigned child_cnt() const override
    {
      return 0;
    }
    virtual unsigned param_cnt() const override
    {
      return 6;
    }
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      
      glm::vec3 size = scene_bbox.max_pos-scene_bbox.min_pos;
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "leg_radius"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "leg_height"});
      params.push_back({5,0.01f*size.z,size.z, ParameterType::DIFFERENTIABLE, "chair_width"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "chair_thickness"});
      params.push_back({5,0.01f*size.x,size.x, ParameterType::DIFFERENTIABLE, "chair_length"});
      params.push_back({5,0.01f*size.y,size.y, ParameterType::DIFFERENTIABLE, "chair_height"});

      return params;
    }
  protected:
    std::vector<std::vector<float>> param_structure() const
    {
      return {{6, ParamNode::CONST, 0}, {7, ParamNode::NEG}, {6},
              {4, ParamNode::PRIMITIVE, 4}, {3, ParamNode::PRIMITIVE, 3}, {2, ParamNode::PRIMITIVE, 2},

              {8, ParamNode::SUB}, {9, ParamNode::NEG}, {6},
              {3}, {5, ParamNode::PRIMITIVE, 5}, {2},

              {10, ParamNode::SUB}, {1, ParamNode::PRIMITIVE, 1}, {11, ParamNode::SUB},
              {1}, {0, ParamNode::PRIMITIVE, 0},

              {12, ParamNode::NEG}, {1}, {11},
              {1}, {0},
              
              {10}, {1}, {13, ParamNode::NEG},
              {1}, {0},
              
              {12}, {1}, {13},
              {1}, {0},
              
              {3},
              {3}, {4},
              {5},
              {4}, {0},
              {2}, {0},
              {10},
              {11}};
    }
    UPGStructure structure() const
    {
      return {{SdfNodeType::OR, SdfNodeType::OR, SdfNodeType::MOVE, SdfNodeType::BOX, SdfNodeType::MOVE, 
               SdfNodeType::BOX, SdfNodeType::OR, SdfNodeType::OR, SdfNodeType::MOVE, SdfNodeType::CYLINDER,
               SdfNodeType::MOVE, SdfNodeType::CYLINDER, SdfNodeType::OR, SdfNodeType::MOVE, SdfNodeType::CYLINDER, 
               SdfNodeType::MOVE, SdfNodeType::CYLINDER}};;
    }
  };
  //test end

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
      SdfNode *node = create_node(n);
      all_nodes.push_back(std::unique_ptr<SdfNode>(node));
      desc.add_parameters(i, node->get_node_name(), node->get_parameters_block(scene_bbox));
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

  std::vector<SdfNodeProperties> node_properties = 
  {
    {SdfNodeType::UNDEFINED  , "UNDEFINED"  , 0, 0, nullptr},
    {SdfNodeType::SPHERE     , "Sphere"     , 1, 0, {[]() -> SdfNode* {return new SphereSdfNode;}}},
    {SdfNodeType::MOVE       , "Move"       , 3, 1, {[]() -> SdfNode* {return new MoveSdfNode;}}},
    {SdfNodeType::OR         , "Or"         , 0, 2, {[]() -> SdfNode* {return new OrSdfNode;}}},
    {SdfNodeType::BOX        , "Box"        , 3, 0, {[]() -> SdfNode* {return new BoxSdfNode;}}},
    {SdfNodeType::CYLINDER   , "Cylinder"   , 2, 0, {[]() -> SdfNode* {return new CylinderSdfNode;}}},
    {SdfNodeType::ROUNDED_BOX, "Rounded Box", 4, 0, {[]() -> SdfNode* {return new RoundBoxSdfNode;}}},
    {SdfNodeType::PRISM      , "Prism"      , 3, 0, {[]() -> SdfNode* {return new PrismSdfNode;}}},
    {SdfNodeType::CONE       , "Cone"       , 4, 0, {[]() -> SdfNode* {return new ConeSdfNode;}}},
    {SdfNodeType::AND        , "And"        , 0, 2, {[]() -> SdfNode* {return new AndSdfNode;}}},
    {SdfNodeType::SUBTRACT   , "Subtract"   , 0, 2, {[]() -> SdfNode* {return new SubtractSdfNode;}}},
    {SdfNodeType::ROTATE     , "Rotate"     , 4, 1, {[]() -> SdfNode* {return new RotateSdfNode;}}},
    {SdfNodeType::SCALE      , "Scale"      , 1, 1, {[]() -> SdfNode* {return new ScaleSdfNode;}}},
    {SdfNodeType::CHAIR      , "Chair"      , VARIABLE_PARAM_COUNT, VARIABLE_CHILD_COUNT, nullptr},
    {SdfNodeType::GRID       , "Grid"       , VARIABLE_PARAM_COUNT, 0, {[]() -> SdfNode* {return new GridSdfNode(32, AABB({-1,-1,-1},{1,1,1}));}}},
    {SdfNodeType::NEURAL     , "Neural"     , VARIABLE_PARAM_COUNT, 0, nullptr},
  };

  const SdfNodeProperties &get_sdf_node_properties(uint16_t type)
  {
    assert(type < SdfNodeType::NODE_TYPES_COUNT);
    assert(node_properties.size() == SdfNodeType::NODE_TYPES_COUNT);
    return node_properties[type];
  }

  const SdfNodeProperties &get_sdf_node_properties(SdfNodeType::Type type)
  {
    assert(node_properties.size() == SdfNodeType::NODE_TYPES_COUNT);
    return node_properties[(int)type];
  }

  SdfNode *create_node(SdfNodeType::Type type)
  {
    assert(node_properties.size() == SdfNodeType::NODE_TYPES_COUNT);
    return node_properties[(int)type].default_constructor();
  }

  SdfNode *create_node(uint16_t type)
  {
    assert(type < SdfNodeType::NODE_TYPES_COUNT);
    assert(node_properties.size() == SdfNodeType::NODE_TYPES_COUNT);
    return node_properties[type].default_constructor();
  }

  int get_position_index(const UPGStructure &structure, int start, int p_start)
  {
    int pi = p_start;
    int si = start;
    while (si < structure.s.size() && get_sdf_node_properties(structure.s[si]).children == 1)
    {
      if (structure.s[si] == SdfNodeType::MOVE)
        return pi;
      si++;
      pi += get_sdf_node_properties(structure.s[si]).param_count;
    }
    return -1;
  }

  int parse_node_rec(const UPGStructure &structure, std::vector<UPGPart> &parts, bool merge_level, int start)
  {

    bool is_merge = merge_level && (structure.s[start] == SdfNodeType::OR);
    int power = get_sdf_node_properties(structure.s[start]).children;
    int pos = start + 1;

    for (int i = 0; i < power; i++)
    {
      int fin_pos = parse_node_rec(structure, parts, is_merge, pos);
      if (is_merge && structure.s[pos] != SdfNodeType::OR)
      {
        int p_st = parts.empty() ? 0 : parts.back().p_range.second;
        int p_cnt = 0;
        for (int p = pos; p < fin_pos; p++)
          p_cnt += get_sdf_node_properties(structure.s[p]).param_count;
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
      cnt += get_sdf_node_properties(node_id).param_count;
    return cnt;
  }

  std::vector<UPGPart> get_sdf_parts(const UPGStructure &structure)
  {
    //this works only if number of parameters and childer for each node is determined by its structure
    for (auto &s : structure.s)
    {
      assert(get_sdf_node_properties(s).param_count != VARIABLE_PARAM_COUNT);
      assert(get_sdf_node_properties(s).children != VARIABLE_CHILD_COUNT);
    }
    std::vector<UPGPart> parts;
    parse_node_rec(structure, parts, true, 0);
    if (parts.empty())
      parts.push_back(UPGPart({0, (int)structure.s.size()}, {0, total_param_count(structure)}, get_position_index(structure, 0, 0)));
    //for (auto &g : parts)
    //  logerr("part [%d %d][%d %d] %d", (int)g.s_range.first, (int)g.s_range.second, g.p_range.first, g.p_range.second, g.position_index);
    return parts;
  }
}