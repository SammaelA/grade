#pragma once
#include <vector>
#include "common_utils/blk.h"
#include "common_utils/utility.h"
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
      float value; //if type is CONST, value is valid
      float min_val;
      float max_val;
      ParameterType type;
      std::string name;
    };
    struct ParamBlock
    {
      std::vector<Param> p;
      std::string name;
    };
    void add_parameters(unsigned node_id, const std::string &node_name, const std::vector<Param> &params)
    {
      block_params[node_id] = {params, node_name};
      prepared = false;
    }
    void remove_parameters(unsigned node_id)
    {
      block_params.erase(node_id);
      prepared = false;
    }
    void add(const ParametersDescription &desc)
    {
      for (const auto &p : desc.block_params)
      {
        if (block_params.find(p.first) == block_params.end())
          block_params[p.first] = p.second;
        else
          logerr("ParameterDescription: trying to merge with description that has the same parameters block %d", p.first);
      }
      prepared = false;
    }
    void print_info()
    {
      int total_params = 0;
      int diff_params = 0;
      int mutable_params = 0;
      int const_params = 0;
      
      for (const auto &b : block_params)
      {
        for (const auto &p : b.second.p)
        {
          total_params++;
          diff_params += (p.type == ParameterType::DIFFERENTIABLE);
          const_params += (p.type == ParameterType::CONST);
          mutable_params += (p.type == ParameterType::MUTABLE_BOOL ||
                             p.type == ParameterType::MUTABLE_INT  ||
                             p.type == ParameterType::MUTABLE_FLOAT);
        }
      }
      debug("=============================\n");
      debug("Parameter description info:\n");
      debug("Total parameters: %d\n", total_params);
      debug("Differentiable parameters: %d\n", diff_params);
      debug("Mutable parameters: %d\n", mutable_params);
      debug("Const parameters: %d\n", const_params);
      for (const auto &b : block_params)
      {
        debug("\tnode %u:\"%s\"", b.first, b.second.name.c_str());
        for (const auto &p : b.second.p)
        {
          switch (p.type)
          {
          case ParameterType::UNKNOWN :
            debug("\t\tu  : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name);
            break;
          case ParameterType::DIFFERENTIABLE :
            debug("\t\td  : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name);
            break;
          case ParameterType::MUTABLE_BOOL :
            debug("\t\tmb : in {true, false}   %s\n", p.name);
            break;
          case ParameterType::MUTABLE_FLOAT :
            debug("\t\tmf : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name);
            break;
          case ParameterType::MUTABLE_INT :
            debug("\t\tmi : in [%8d, %8d)   %s\n", (int)p.min_val, (int)p.max_val, p.name);
            break;
          case ParameterType::CONST :
            debug("\t\tc : %8.4f   %s\n", p.value, p.name);
            break;
          default:
            break;
          }
        }
      }
      debug("=============================\n");
    }
  private:
    bool prepared = false;
    std::unordered_map<unsigned, ParamBlock> block_params;
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