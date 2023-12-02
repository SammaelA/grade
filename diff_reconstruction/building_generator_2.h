#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  template <typename real>
  PartOffsets create_building_2t(const std::vector<real> &params, std::vector<real> &model, ModelQuality quality);

  PartOffsets create_building_2(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality);
  PartOffsets create_building_2f(const std::vector<float> &params, std::vector<float> &model, ModelQuality quality);
}