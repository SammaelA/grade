#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  void create_simple_table(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality);
}