#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include <memory>

namespace upg
{
  class SdfNode
  {
  protected:
    unsigned ID;
    std::string name;
    std::span<const float> p;
  public:
    SdfNode(unsigned id) { ID = id; }
    virtual ~SdfNode() = default;
    unsigned get_ID() const { return ID; }
    std::string get_node_name() const { return name; }
    void set_param_span(std::span<my_float> s) { p = s; }

    virtual float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const = 0;
    virtual bool add_child(SdfNode *node) = 0;// returns the availability of free space
    virtual unsigned param_cnt() const = 0;
    virtual unsigned child_cnt() const = 0;
    virtual std::vector<const SdfNode *> get_children() const = 0;
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const = 0;
    virtual AABB get_bbox() const = 0;
  };

  class ProceduralSdf
  {
  public:
    ProceduralSdf() = default;
    ProceduralSdf(const SdfNode *_root, const UPGStructure &_structure,
                std::span<const std::unique_ptr<SdfNode>> _all_nodes):
      all_nodes(_all_nodes)
    {
      root = _root;
      structure = _structure;
    }
    float get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                               std::vector<float> *ddist_dpos = nullptr) const
    {
      return root->get_distance(pos, ddist_dp, ddist_dpos);
    }
    const SdfNode *root = nullptr;
    const std::span<const std::unique_ptr<SdfNode>> all_nodes;
    UPGStructure structure;
  };

  class SdfGenInstance : public UniversalGenInstance
  {
  public:
    //to better perform optimization on different scales, we should better
    //determine the size of scene we are working on.
    //BBox affects only ParametersDescription, not the generation itself.
    static void set_scene_bbox(AABB bbox) {scene_bbox = bbox;}

    SdfGenInstance(const UPGStructure &structure);
    virtual void recreate(const UPGStructure &structure) override;
    ProceduralSdf generate(std::span<const float> parameters);
    ParametersDescription desc;

  private:
    std::vector<std::unique_ptr<SdfNode>> all_nodes;
    SdfNode *root;
    //spans from inputParams points to this container
    //put raw parameters list here to generate
    //DO NOT change size of this vector
    std::vector<my_float> all_params;
    UPGStructure structure;

    static AABB scene_bbox;
  };

  SdfNode *sdf_node_by_node_type_id(uint16_t num, unsigned id);
  std::vector<UPGPart> get_sdf_parts(const UPGStructure &structure);
}