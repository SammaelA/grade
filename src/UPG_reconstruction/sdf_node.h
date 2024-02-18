#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include <memory>

namespace upg
{
  class ProceduralSdf;
  class SdfNode
  {
  public:
    enum NodeType
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
      GRID,
      NEURAL,
      NODE_TYPES_COUNT
    };

    friend class ProceduralSdf;

    SdfNode(unsigned id) { ID = id; }
    virtual ~SdfNode() = default;
    unsigned get_ID() const { return ID; }
    std::string get_node_name() const { return name; }
    void set_param_span(std::span<my_float> s, unsigned offset) 
    { 
      p = s; 
      p_offset = offset;
    }

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,     //3*batch_size: p0.x, p0.y, p0.z, p1.x, ...
                                    float *      distances,     //batch_size
                                    float *      ddist_dparams, //Nparams*batch_size, one array for all nodes
                                    float *      ddist_dpos,    //3*batch_size: ddist_p0.x, ddist_p0.y, ddist_p0.z, ddist_p1.x, ...
                            std::vector<float> & stack,         //resizable chunk of memory for temporary data (should only grow in size)
                                    unsigned     stack_head) const   //current position of stack head (0 initially, grows up)
    {
      logerr("NOT IMPLEMENTED");
    }
    virtual bool add_child(SdfNode *node) = 0;// returns the availability of free space
    virtual unsigned param_cnt() const = 0;
    virtual unsigned child_cnt() const = 0;
    virtual std::vector<const SdfNode *> get_children() const = 0;
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const = 0;
    virtual AABB get_bbox() const = 0;
  //protected:
    unsigned ID;
    unsigned p_offset = 0;
    unsigned subgraph_param_cnt = 0;
    std::vector<SdfNode *> subgraph;
    std::string name;
    std::span<const float> p;
  };

  class PrimitiveSdfNode : public SdfNode
  {
  public:
    PrimitiveSdfNode(unsigned id) : SdfNode(id) {}
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

  SdfNode *sdf_node_by_node_type_id(uint16_t num, unsigned id);
  std::vector<UPGPart> get_sdf_parts(const UPGStructure &structure);
}