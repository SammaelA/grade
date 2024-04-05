#include "sdf_octree.h"

namespace upg
{
  void OctreeSdfNode::get_distance_batch(unsigned     batch_size,
                                       float *const positions,    
                                       float *      distances,
                                       float *      ddist_dparams,
                                       float *      ddist_dpos,
                               std::vector<float> & stack,
                                       unsigned     stack_head) const
  {
    if (ddist_dparams)
    {
      logerr("Sdf OCtree differentiation is not implemented!");
      assert(false);
    }
    else
    {
      for (int i=0;i<batch_size;i++)
        distances[i] = octree.sample(float3(positions[3*i+0], positions[3*i+1], positions[3*i+2]));
    }
  }

  std::vector<ParametersDescription::Param> OctreeSdfNode::get_parameters_block(AABB scene_bbox) const
  {
    std::vector<ParametersDescription::Param> params;
    params.push_back({1,-1.0f,1.0f, ParameterType::ARRAY, "data", nodes_limit});
    return params;
  }

  void OctreeSdfNode::construct(std::function<float(const float3 &)> sdf)
  {
    SparseOctreeSettings settings;
    settings.depth = 8;
    BlockSparseOctree<float> bso;
    octree.construct_bottom_up(sdf, settings);
  }
}