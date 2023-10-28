#pragma once
#include <vector>
#include "common_utils/blk.h"
#include "common_utils/utility.h"
#include "third_party/span.h"
#include <unordered_map>

namespace upg
{
  struct UPGStructure
  {

  };
  struct UPGReconstructionResult
  {
    float quality = 0;
    UPGStructure structure;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk);
};