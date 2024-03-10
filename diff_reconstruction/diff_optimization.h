#pragma once
#include <vector>

#define USE_OPEN_CV 0

struct Block;
class MitsubaInterface;
namespace dopt
{
  struct OptimizationResult
  {
    std::vector<float> best_params;
    float best_err;
    int total_iters;
  };
  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi);
}