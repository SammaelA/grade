#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  //create model of a building with given parameters (buildings procedural generator)
  PartOffsets create_building(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality);
}