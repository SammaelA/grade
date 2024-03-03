#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include <memory>
#include <functional>

namespace upg
{
  class SdfNode;
  struct SdfNodeType
  {
    enum Type
    {
      UNDEFINED,
      SPHERE,
      MOVE,
      OR,
      BOX,
      CYLINDER,
      ROUNDED_BOX,
      PRISM,
      CONE,
      AND,
      SUBTRACT,
      ROTATE,
      SCALE,
      CHAIR,
      CROTATE,
      ROUND,
      CIRCLE,
      QUAD,
      EXTRUSION,
      SCALE2D,
      MOVE2D,
      ROTATE2D,
      CROTATE2D,
      POLY8,
      GRID_16,
      GRID_32,
      GRID_64,
      GRID_128,
      GRID_256,
      NEURAL_TINY,
      NEURAL_SMALL,
      NEURAL_MEDIUM,
      NEURAL_LARGE,
      NEURAL_HUGE,
      NODE_TYPES_COUNT
    };
  };

  enum class SdfNodeClass
  {
    PRIMITIVE,
    TRANSFORM,
    COMBINE,
    COMPLEX,
    GRID,
    NEURAL,
    TRANSFORM2D,
    PRIMITIVE2D,
    OTHER
  };

  static constexpr int VARIABLE_PARAM_COUNT = -1;
  static constexpr int VARIABLE_CHILD_COUNT = -1;

  struct SdfNodeProperties
  {
    SdfNodeType::Type type;
    std::string name;
    SdfNodeClass node_class;
    int param_count;
    int children;
    std::function<SdfNode *(SdfNodeType::Type)> default_constructor;
  };

  const SdfNodeProperties &get_sdf_node_properties(uint16_t id);
  const SdfNodeProperties &get_sdf_node_properties(SdfNodeType::Type type);
  
  SdfNode *create_node(uint16_t id);
  SdfNode *create_node(SdfNodeType::Type type);

  class ProceduralSdf;
  class SdfNode
  {
  public:
    friend class ProceduralSdf;

    SdfNode(const SdfNodeType::Type &_type):
      type(_type)
    {

    }
    virtual ~SdfNode() = default;
    void set_param_span(std::span<float> s, unsigned offset) 
    { 
      p = s; 
      p_offset = offset;
    }
    virtual void set_subgraph_params_cnt_rec() const
    {
      unsigned cnt = param_cnt();
      for (auto *c : get_children())
      {
        c->set_subgraph_params_cnt_rec();
        subgraph.insert(subgraph.end(), c->subgraph.begin(), c->subgraph.end());
        cnt += c->subgraph_param_cnt;
      }
      subgraph_param_cnt = cnt;
      subgraph.push_back(this);
    }

    //method is called every time when node's parameters 
    //are changed. Some nodes overrides it to pre-calculate 
    //some values that depend on parameters
    virtual void on_params_change() const
    {
      for (auto i : get_children())
      {
        i->on_params_change();
      }
    }
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,     //3*batch_size: p0.x, p0.y, p0.z, p1.x, ...
                                    float *      distances,     //batch_size
                                    float *      ddist_dparams, //Nparams*batch_size, one array for all nodes
                                    float *      ddist_dpos,    //3*batch_size: ddist_p0.x, ddist_p0.y, ddist_p0.z, ddist_p1.x, ...
                            std::vector<float> & stack,         //resizable chunk of memory for temporary data (should only grow in size)
                                    unsigned     stack_head) const = 0;//current position of stack head (0 initially, grows up)

    virtual bool add_child(SdfNode *node) = 0;// returns the availability of free space
    virtual unsigned param_cnt() const = 0;
    virtual unsigned child_cnt() const = 0;
    virtual std::vector<const SdfNode *> get_children() const = 0;
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const = 0;
    virtual AABB get_bbox() const = 0;
    SdfNodeType::Type get_type() const { return type; }
  //protected:
    unsigned p_offset = 0;
    mutable unsigned subgraph_param_cnt = 0;
    mutable std::vector<const SdfNode *> subgraph;
    std::span<const float> p;
    SdfNodeType::Type type = SdfNodeType::UNDEFINED;
  };

  class PrimitiveSdfNode : public SdfNode
  {
  public:
    PrimitiveSdfNode(const SdfNodeType::Type &_type) : SdfNode(_type) {}
    virtual bool add_child(SdfNode *node) override { return false; }
    virtual unsigned child_cnt() const override { return 0; }
    virtual std::vector<const SdfNode *> get_children() const override { return {}; }
  };

  class ProceduralSdf : public UniversalGenInstance
  {
  public:
    //to better perform optimization on different scales, we should better
    //determine the size of scene we are working on.
    //BBox affects only ParametersDescription, not the generation itself.
    static void set_scene_bbox(AABB bbox) {scene_bbox = bbox;}

    ProceduralSdf(const UPGStructure &structure);
    ProceduralSdf(const ProceduralSdf &sdf);
    virtual void recreate(const UPGStructure &structure) override;
    void set_parameters(std::span<const float> parameters);
    ParametersDescription desc;

    float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const;
    void get_distance_batch(unsigned     batch_size,
                            float *const positions,     //3*batch_size: p0.x, p0.y, p0.z, p1.x, ...
                            float *      distances,     //batch_size
                            float *      ddist_dparams, //Nparams*batch_size, one array for all nodes
                            float *      ddist_dpos) const;//3*batch_size: ddist_p0.x, ddist_p0.y, ddist_p0.z, ddist_p1.x, ...
    AABB get_bbox() const
    {
      return root->get_bbox();
    }
  //private:
    std::vector<std::unique_ptr<SdfNode>> all_nodes;
    SdfNode *root;
    //spans from inputParams points to this container
    //put raw parameters list here to generate
    //DO NOT change size of this vector
    std::vector<float> all_params;
    mutable std::vector<float> stack;
    mutable std::vector<float> ddist_dparams_transp;
    UPGStructure structure;

    static AABB scene_bbox;
    
  };

  std::vector<UPGPart> get_sdf_parts(const UPGStructure &structure);
  bool is_struct_correct(const UPGStructure &structure);
}