#pragma once
#include "differentiable_generators.h"

namespace dgen
{
  double iou3d(const std::vector<dfloat> &model1, const std::vector<dfloat> &model2);
}