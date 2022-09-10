#pragma once
#include <vector>
#include "vectors.h"
struct ComplexModel;
struct Block;
namespace dgen
{
  void dgen_test(std::vector<float> &model);
  bool create_model_from_block(Block &bl, ComplexModel &mod);
};