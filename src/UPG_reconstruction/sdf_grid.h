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
    GridSdfNode(const SdfNodeType::Type &_type, unsigned _grid_size) : PrimitiveSdfNode(_type)
    { 
      grid_size = _grid_size;
      grid_size_f = grid_size;
    }
    virtual ~GridSdfNode() = default;

    float get_distance(const float3 &pos, std::span<float> ddist_dp, 
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
    //grid SDF is always set in a unit cube and transformed by other nodes if needed
    virtual AABB get_bbox() const override { return AABB({-1,-1,-1},{1,1,1}); }
    void set_voxel(const float3 &pos, float distance);
    void set_voxel(const uint3 &vox, float distance);

    void save(const std::string &path) const;
    void load(const std::string &path);
  protected:
    float sample_bilinear(const float3 &pos, std::vector<float> *ddist_dp = nullptr, 
                          std::vector<float> *ddist_dpos = nullptr) const;

    unsigned grid_size;
    float grid_size_f;
    float3 bbox_size = float3(2,2,2);
  };
}