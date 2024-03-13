#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include "common_utils/sparse_octree.h"
#include "sdf_node.h"
#include <memory>
#include <functional>

namespace upg
{
  class OctreeSdfNode : public PrimitiveSdfNode
  {
  public:
    OctreeSdfNode(const SdfNodeType::Type &_type, unsigned _nodes_limit) : PrimitiveSdfNode(_type)
    { 
      nodes_limit = _nodes_limit;
    }
    virtual ~OctreeSdfNode() = default;

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,    
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override;

    virtual unsigned param_cnt() const override { return nodes_count; } 
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override;
    //grid SDF is always set in a unit cube and transformed by other nodes if needed
    virtual AABB get_bbox() const override { return AABB({-1,-1,-1},{1,1,1}); }
    void construct(std::function<float(const float3 &)> sdf);
  //protected:
    float sample_bilinear(const float3 &pos, std::vector<float> *ddist_dp = nullptr, 
                          std::vector<float> *ddist_dpos = nullptr) const;

    unsigned nodes_limit = 1;
    unsigned nodes_count = 1;
    float3 bbox_size = float3(2,2,2);
    SparseOctree octree;
  };
}