#pragma once
#include <vector>
#include <map>
#include "p_v_structs.h"
#include "common_utils/utility.h"
namespace upg
{
  struct UniversalGenMesh
  {
    //triangle mesh pos.size()%9 == 0
    //norm and tc can be empty
    std::vector<float> pos; //vec3
    std::vector<float> norm; //vec3
    std::vector<float> tc; //vec2
  };

  void add_tri_data(upg::vec3 point, upg::vec3 n, upg::vec2 tex, UniversalGenMesh &mesh);
  void add_point_data(upg::vec3 point, UniversalGenMesh &mesh);

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
    int get_total_params_count() const
    {
      return total_params_count;
    }
  private:
    int total_params_count = 0;
    std::map<unsigned, ParamBlock> block_params;
  };
}