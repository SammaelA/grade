#pragma once
#include "reconstruction.h"
#include "tinyEngine/texture.h"

namespace upg
{
  //structure that contains all data for one given view, like mask and maybe some 
  //depth information
  struct ReferenceView
  {
    Texture mask;
  };
  class UniversalGen
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
        debug("  node %u: \"%s\"\n", b.first, b.second.name.c_str());
        for (const auto &p : b.second.p)
        {
          switch (p.type)
          {
          case ParameterType::UNKNOWN :
            debug("    u  : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::DIFFERENTIABLE :
            debug("    d  : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::MUTABLE_BOOL :
            debug("    mb : in {true, false}   %s\n", p.name.c_str());
            break;
          case ParameterType::MUTABLE_FLOAT :
            debug("    mf : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::MUTABLE_INT :
            debug("    mi : in [%8d, %8d)   %s\n", (int)p.min_val, (int)p.max_val, p.name.c_str());
            break;
          case ParameterType::CONST :
            debug("    c : %8.4f   %s\n", p.value, p.name.c_str());
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
}