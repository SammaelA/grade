#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  //create model of a cup with given parameters (dishes procedural generator)
  void create_cup(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality);
  //estimates how weird are parametes of this cup (0 for good parameters)
  dfloat parameters_cup_reg(const std::vector<dfloat> &params);
}