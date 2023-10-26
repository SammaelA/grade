#pragma once
#include <vector>
#include "common_utils/blk.h"
#include "third_party/span.h"
#include <unordered_map>

namespace upg
{
  class UniversalGen
  {

  };
  struct UPGStructure
  {

  };
  enum class ParameterType
  {
    UNKNOWN,
    DIFFERENTIABLE,
    MUTABLE_FLOAT,
    MUTABLE_INT,
    MUTABLE_BOOL,
    CONST
  };
  struct ParametersDescription
  {
    struct Param
    {
      float min_val;
      float max_val;
      ParameterType type;
      std::string name;
    };
    void add_parameters(unsigned node_id, const std::vector<Param> &params)
    {
      param_by_block[node_id] = params;
      prepared = false;
    }
    void remove_parameters(unsigned node_id)
    {
      param_by_block.erase(node_id);
      prepared = false;
    }
  private:
    bool prepared = false;
    std::unordered_map<unsigned, std::vector<Param>> param_by_block;
  };
  
  struct UPGNodeInputParameters
  {

    std::span<const float> diff_params;
    std::span<const float> not_diff_params;
  };
  struct UPGInputParameters
  {
    std::vector<float> values;
    std::unordered_map<unsigned, UPGNodeInputParameters> nodes; //node parameters by id
  };

  struct UPGReconstructionResult
  {
    float quality = 0;
    UPGStructure structure;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk);
};