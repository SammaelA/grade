#pragma once
#include <vector>
#include "common_utils/blk.h"

namespace upg
{
  class UniversalGen
  {

  };
  struct UPGStructure
  {

  };
  struct UPGParameters
  {

  };
  struct UPGCameraParameters
  {
    
  };
  struct UPGReconstructionResult
  {
    float quality = 0;
    UPGStructure structure;
    UPGParameters parameters;
    UPGCameraParameters camera_parameters;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk);
};