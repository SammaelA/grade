#pragma once
#include "differentiable_generators.h"
struct ComplexModel;
struct Block;
namespace dgen
{
  void print_jackobian(const std::vector<float> &jac, int x_n, int y_n, int lines = 100);
  void print_model(const std::vector<float> &res);
  void dgen_test(std::string generator_name, std::vector<float> &params, std::vector<float> &model);
  void dgen_test_internal(std::vector<float> &model, generator_func func, const std::vector<float> &check_params, 
                          const std::vector<float> &params, std::vector<float> *jacobian = nullptr);
  bool check_stability(generator_func func, const std::vector<float> &params, int iterations);
  bool check_robustness(generator_func func, const std::vector<float> &params_min, const std::vector<float> &params_max, int iterations);
  bool create_model_from_block(Block &bl, ComplexModel &mod);
  void transform(std::vector<dfloat> &vert, dmat43 mat, int floats_per_vertex = FLOAT_PER_VERTEX, int pos_start = 0, int norm_start = 3);
  void transform(std::vector<dvec3> &verts, dmat43 mat);
  void transform_by_scene_parameters(std::vector<dgen::dfloat> &params, int offset, std::vector<dgen::dfloat> &model);
  dfloat smoothmax(dfloat a, dfloat b, float alpha = 16);
  dfloat smoothmin(dfloat a, dfloat b, float alpha = 16);
  dfloat smoothclamp(dfloat a, dfloat from, dfloat to, float alpha = 16);
  void add_vertex(std::vector<dfloat> &vert, int n, const dvec3 &pos, const dvec3 &norm, const dvec2 &tc, bool only_pos);
  //returns 0 is each parameters is in [min+edge_size, max-edge_size] interval
  dfloat parameters_limits_reg(const std::vector<dfloat> &params, const std::vector<float> &params_min, const std::vector<float> &params_max,
                               float edge_size = 0.01);
  dfloat default_parameters_reg(const std::vector<dfloat> &params);
  dfloat default_model_reg(const std::vector<dfloat> &params, const std::vector<dfloat> &model);
};