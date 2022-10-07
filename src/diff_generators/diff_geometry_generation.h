#pragma once
#include <vector>
#include "vectors.h"
struct ComplexModel;
struct Block;
namespace dgen
{
  #define FLOAT_PER_VERTEX (3+3+2) //vec3 pos, vec3 norm, vec2 tc
  void print_jackobian(const std::vector<float> &jac, int x_n, int y_n, int lines = 100);
  void print_model(const std::vector<float> &res);
  void dgen_test(std::vector<float> &model);
  void dgen_test_internal(std::vector<float> &model, const std::vector<float> &check_params, const std::vector<float> &params,
                          std::vector<float> *jacobian = nullptr);
  bool check_stability(const std::vector<float> &params, int iterations);
  bool check_robustness(const std::vector<float> &params_min, const std::vector<float> &params_max, int iterations);
  bool create_model_from_block(Block &bl, ComplexModel &mod);
  void transform(std::vector<dfloat> &vert, dmat43 mat, int floats_per_vertex = FLOAT_PER_VERTEX, int pos_start = 0, int norm_start = 3);
  void create_cup(std::vector<dfloat> &params, std::vector<dfloat> &model);
  void create_plate(std::vector<dfloat> &params, std::vector<dfloat> &model);
};