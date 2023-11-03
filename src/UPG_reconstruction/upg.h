#pragma once
#include <vector>
#include "common_utils/blk.h"
#include "common_utils/utility.h"
#include "third_party/span.h"
#include "tinyEngine/model.h"
#include <map>

namespace upg
{
  struct UPGStructure
  {
    std::vector<uint16_t> s;
  };
  struct UPGParametersRaw
  {
    std::vector<float> p;
  };
  struct UPGReconstructionResult
  {
    float quality = 0;
    UPGStructure structure;
    UPGParametersRaw parameters;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk);
  bool create_model_from_block(const Block &bl, ComplexModel &out_mod);
  bool create_model(const UPGStructure &structure, const UPGParametersRaw &params,
                    ComplexModel &out_mod);
};