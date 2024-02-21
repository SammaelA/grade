#pragma once
#include "upg.h"
#include "generation_common.h"
#include "common_utils/bbox.h"
#include "sdf_node.h"
#include "neuralCore/siren.h"
#include <memory>

namespace upg
{
  class NeuralSdfNode : public PrimitiveSdfNode
  {
  public:
    NeuralSdfNode(unsigned hidden_layers, unsigned layer_size, const AABB &region_bbox) : 
    PrimitiveSdfNode(),
    network(nn::Siren::Type::SDF, hidden_layers, layer_size)
    { 
      network.initialize();
      network.set_batch_size_for_evaluate(base_batch_size);
      params_count = network.params_count();
      bbox = region_bbox.expand(1/1.1f);
      name = "Neural";
    }
    virtual ~NeuralSdfNode() = default;

    virtual void get_distance_batch(unsigned     batch_size,
                                    float *const positions,    
                                    float *      distances,
                                    float *      ddist_dparams,
                                    float *      ddist_dpos,
                            std::vector<float> & stack,
                                    unsigned     stack_head) const override
    {
      if (!weights_set)
      {
        weights_set = true;
        network.initialize_with_weights(p.data());
      }
      if (ddist_dparams == nullptr && ddist_dpos == nullptr)
        network.evaluate(positions, distances, batch_size);
      else
        logerr("NeuralSdfNode: differentiation of neural SDF is not implemented!");
    }

    virtual unsigned param_cnt() const override { return params_count; } 
    virtual std::vector<ParametersDescription::Param> get_parameters_block(AABB scene_bbox) const override
    {
      std::vector<ParametersDescription::Param> params;
      params.push_back({1,-1.0f,1.0f, ParameterType::ARRAY, "data", params_count});
      return params;
    }
    virtual AABB get_bbox() const override { return bbox; };
  protected:

    unsigned params_count;
    unsigned base_batch_size = 1024;
    mutable nn::Siren network;
    mutable bool weights_set = false;
    AABB bbox;
  };
}