#pragma once
#include "reconstruction_sdf.h"

namespace upg
{
  std::vector<UPGReconstructionResult> reconstruction_graph_based(Block *step_blk, const std::vector<glm::vec3> &points,
                                                                  const std::vector<float> &distances);
}