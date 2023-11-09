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
    //loss of optimization function that was achieved on this model 
    float loss_optimizer = 0;
    //similarity function (default sil-PSNR) between reference images and reconstructed model rendered from the same camera positions
    float quality_ir = 0;
    //similarity function (default sil-PSNR on 64 views) between reference model and reconstructed model. Only for synthetic reference
    float quality_synt = 0;
    UPGStructure structure;
    UPGParametersRaw parameters;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk);
  bool create_model_from_block(const Block &bl, ComplexModel &out_mod);
  bool create_model(const UPGStructure &structure, const UPGParametersRaw &params,
                    ComplexModel &out_mod);

  void perform_tests(const Block &blk);
};