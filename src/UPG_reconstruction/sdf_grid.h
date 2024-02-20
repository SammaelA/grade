#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include "sdf_node.h"
#include <memory>

namespace upg
{
  class GridSdfNode : public PrimitiveSdfNode
  {
  public:
    GridSdfNode(unsigned _grid_size, const AABB &region_bbox) : PrimitiveSdfNode()
    { 
      grid_size = _grid_size;
      bbox = region_bbox;
      bbox_size = bbox.size();
      grid_size_f = grid_size;
      name = "Grid";
    }
    virtual ~GridSdfNode() = default;

    float get_distance(const glm::vec3 &pos, std::span<float> ddist_dp, 
                       std::span<float> ddist_dpos) const;
    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,    
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override;

    virtual unsigned param_cnt() const override { return grid_size*grid_size*grid_size; } 
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override;
    virtual AABB get_bbox() const override { return bbox; };
    void set_voxel(const glm::vec3 &pos, float distance);
    void set_voxel(const glm::uvec3 &vox, float distance);
  protected:
    float sample_bilinear(const glm::vec3 &pos, std::vector<float> *ddist_dp = nullptr, 
                          std::vector<float> *ddist_dpos = nullptr) const;

    unsigned grid_size;
    float grid_size_f;
    AABB bbox;
    glm::vec3 bbox_size;
  };
}