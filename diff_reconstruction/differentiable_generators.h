#pragma once
#include <vector>
#include "cpp_ad_vectors.h"
#include <functional>
#include "common_utils/utility.h"

namespace dgen
{
  #define FLOAT_PER_VERTEX (3+3+2) //vec3 pos, vec3 norm, vec2 tc
  struct ModelLayout
  {
    //default layout is (pos.x, pos.y, pos.z, norm.x, norm.y, norm.z, tc.x, tc.y)
    //if some value is -1, it means that model does not have such component in vertex
    ModelLayout(): ModelLayout(0, 3, 6, 8, 8) {}
    ModelLayout(int _p, int _n, int _tc, int _end, int _offset)
    {
      pos = _p;
      norm = _n;
      tc = _tc;
      end = _end;
      f_per_vert = _offset;
    }
    union
    {
      std::array<int, 4> offsets;
      struct
      {
        int pos;
        int norm;
        int tc;
        int end;
      };
    };
    int f_per_vert = 8;
  };

  struct ModelQuality
  {
    bool create_only_position = false;
    enum Level {LOW, MEDIUM, HIGH, ULTRA} quality_level = MEDIUM;
    
    ModelQuality(bool only_pos, int quality)
    {
      create_only_position = only_pos;
      quality_level = (Level)CLAMP(quality,0,3);
    }
    ModelQuality() : ModelQuality(false, 1) {}
  };

  //composite model is made from several parts
  //this structure is needed to determine what parts
  //the model has and offsets in model data for each part.
  typedef std::vector<std::pair<std::string, int>> PartOffsets;

  PartOffsets simple_mesh();
  
  typedef std::pair<std::vector<float>, PartOffsets> DFModel;

  typedef std::function<PartOffsets(const std::vector<dfloat> &params, std::vector<dfloat> &out_model, ModelQuality)> generator_func;
  typedef std::function<dfloat(const std::vector<dfloat> &params)> params_regularizer_func;
  typedef std::function<dfloat(const std::vector<dfloat> &params, const std::vector<dfloat> &model)> model_regularizer_func;

  struct GeneratorDescription
  {
    std::string name;
    generator_func generator;
    std::function<PartOffsets(const std::vector<float> &params, std::vector<float> &out_model, ModelQuality)> gen_not_diff;
    params_regularizer_func params_regularizer;
    model_regularizer_func model_regularizer;
    std::string generator_description_blk_path;
    std::string presets_blk_path;
  };
  GeneratorDescription get_generator_by_name(std::string name);
};