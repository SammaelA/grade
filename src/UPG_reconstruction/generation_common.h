#pragma once
#include <vector>
#include <map>
#include "common_utils/utility.h"
#include "common_utils/template_vectors.h"
#include "upg.h"
namespace upg
{
  typedef float my_float;

  typedef dgen::g_vec2<my_float> vec2;
  typedef dgen::g_vec3<my_float> vec3;
  typedef dgen::g_vec4<my_float> vec4;
  typedef dgen::g_mat43<my_float> mat43;

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
      float value; //if type is CONST, value is valid, otherwise not.
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
      total_params_count += params.size();
      block_params[node_id] = {params, node_name};
    }
    void remove_parameters(unsigned node_id)
    {
      auto it = block_params.find(node_id);
      if (it != block_params.end())
      {
        total_params_count -= it->second.p.size();
        block_params.erase(node_id);
      }
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
      total_params_count += desc.total_params_count;
    }
    void print_info() const
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
            debug("    u : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::DIFFERENTIABLE :
            debug("    d : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::MUTABLE_BOOL :
            debug("   mb : in {true, false}   %s\n", p.name.c_str());
            break;
          case ParameterType::MUTABLE_FLOAT :
            debug("   mf : in [%8.4f, %8.4f)   %s\n", p.min_val, p.max_val, p.name.c_str());
            break;
          case ParameterType::MUTABLE_INT :
            debug("   mi : in [%8d, %8d)   %s\n", (int)p.min_val, (int)p.max_val, p.name.c_str());
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
    const std::map<unsigned, ParamBlock> &get_block_params() const
    {
      return block_params;
    }
    std::map<unsigned, ParamBlock> &get_block_params()
    {
      return block_params;
    }
    int get_total_params_count() const
    {
      return total_params_count;
    }
  private:
    int total_params_count = 0;
    std::map<unsigned, ParamBlock> block_params;
  };

  class UniversalGenInstance
  {
  public:
    UniversalGenInstance() = default;
    virtual ~UniversalGenInstance() {};
    UniversalGenInstance(const UniversalGenInstance &) = delete;
    UniversalGenInstance &operator=(const UniversalGenInstance&) = delete;
    virtual void recreate(const UPGStructure &structure) = 0;
  };
  
  struct UPGPart
  {
    UPGPart() = default;
    UPGPart(std::pair<int, int> _s, std::pair<int, int> _p, int _position_index)
    {
      s_range = _s;
      p_range = _p;
      position_index = _position_index;
    }
    std::pair<int, int> s_range; // range of indices in structure that corresponds to this part. This range is guaranteed to be a valid structure.
    std::pair<int, int> p_range; // range of indices in parameters array.
    int position_index = -1;     // index of position parameter (first of 3) in parameters array
  };
}