#pragma once
#include "differentiable_generators.h"
struct ComplexModel;
struct Block;
namespace dgen
{
  void print_jackobian(const std::vector<float> &jac, int x_n, int y_n, int lines = 100);
  void print_model(const std::vector<float> &res);
  void dgen_test(std::string generator_name, std::vector<float> &params, dgen::DFModel &model, bool transform_by_scene = false,
                 ModelQuality mq = ModelQuality(false, 1));
  void dgen_test_internal(dgen::DFModel &model, generator_func func, const std::vector<float> &check_params, 
                          const std::vector<float> &params, std::vector<float> *jacobian, bool transform_by_scene,
                          ModelQuality mq);
  bool check_stability(generator_func func, const std::vector<float> &params, int iterations);
  bool check_robustness(generator_func func, const std::vector<float> &params_min, const std::vector<float> &params_max, int iterations);
  bool create_model_from_block(Block &bl, ComplexModel &mod);
  void transform_by_scene_parameters(std::vector<dgen::dfloat> &params, int offset, std::vector<dgen::dfloat> &model);
  void transform_by_scene_parameters(const std::vector<float> &scene_params, std::vector<float> &model);
  dfloat smoothmax(dfloat a, dfloat b, float alpha = 16);
  dfloat smoothmin(dfloat a, dfloat b, float alpha = 16);
  dfloat smoothclamp(dfloat a, dfloat from, dfloat to, float alpha = 16);
  template<typename float_type>
  void add_vertex(std::vector<float_type> &vert, int n, const g_vec3<float_type> &pos, const g_vec3<float_type> &norm, const g_vec2<float_type> &tc, bool only_pos)
  {
    int sz = n*FLOAT_PER_VERTEX;
    vert.resize(sz + FLOAT_PER_VERTEX);
    vert[sz+0] = pos[0];
    vert[sz+1] = pos[1];
    vert[sz+2] = pos[2];
    //if (!only_pos)
    //{
      vert[sz+3] = norm[0];
      vert[sz+4] = norm[1];
      vert[sz+5] = norm[2]; 
      vert[sz+6] = tc[0];
      vert[sz+7] = tc[1];
  }
  //returns 0 is each parameters is in [min+edge_size, max-edge_size] interval
  dfloat parameters_limits_reg(const std::vector<dfloat> &params, const std::vector<float> &params_min, const std::vector<float> &params_max,
                               float edge_size = 0.01);
  dfloat default_parameters_reg(const std::vector<dfloat> &params);
  dfloat default_model_reg(const std::vector<dfloat> &params, const std::vector<dfloat> &model);
  dfloat d_max(dfloat a, dfloat b);
  dfloat d_min(dfloat a, dfloat b);

  template <typename T>
  void transform(std::vector<T> &vert, g_mat43<T> mat, int floats_per_vertex = FLOAT_PER_VERTEX, int pos_start = 0, int norm_start = 3)
  {
    g_mat43<T> norm_mat = transposedInverse3x3(mat);

    if (norm_start >= 0)
    {
      for (int i=0;i<vert.size()/floats_per_vertex;i++)
      {
        mulp(mat, vert[pos_start+floats_per_vertex*i], vert[pos_start+floats_per_vertex*i+1], vert[pos_start+floats_per_vertex*i+2]);
        mulv(norm_mat, vert[norm_start+floats_per_vertex*i], vert[norm_start+floats_per_vertex*i+1], vert[norm_start+floats_per_vertex*i+2]);

        T a = vert[norm_start+floats_per_vertex*i];
        T b = vert[norm_start+floats_per_vertex*i+1];
        T c = vert[norm_start+floats_per_vertex*i+2];
        T len = sqrt(a*a + b*b + c*c) + 1e-18;
        vert[norm_start+floats_per_vertex*i] = a/len;
        vert[norm_start+floats_per_vertex*i+1] = b/len;
        vert[norm_start+floats_per_vertex*i+2] = c/len;
      }
    }
    else
    {
      for (int i=0;i<vert.size()/floats_per_vertex;i++)
      {
        mulp(mat, vert[pos_start+floats_per_vertex*i], vert[pos_start+floats_per_vertex*i+1], vert[pos_start+floats_per_vertex*i+2]);
      }
    }
  }

  template <typename T>
  void transform(std::vector<g_vec3<T>> &verts, g_mat43<T> mat)
  {
    for (auto &vert : verts)
    {
      vert = mulp(mat, vert); 
    }
  }
};