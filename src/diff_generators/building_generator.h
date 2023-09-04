#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  //create model of a building with given parameters (buildings procedural generator)
  template<typename float_type>
  extern PartOffsets create_building(const std::vector<float_type> &params, std::vector<float_type> &model, ModelQuality quality);
}