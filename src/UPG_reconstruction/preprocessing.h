#pragma once
#include "reconstruction_impl.h"

namespace upg
{
  std::vector<ReferenceView> get_reference(const Block &input_blk);
  ReferenceView preprocess_get_reference_view(const Block &view_blk);
  ReferenceView preprocess_get_reference_view_synthetic(const Block &synthetic_reference_blk, const Block &view_blk);
}