#include "diff_optimization.h"
#include "diff_geometry_generation.h"
#include "mitsuba_python_interaction.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"

namespace dopt
{
  void test()
  {
    MitsubaInterface mi;
    mi.init("/home/sammael/grade/scripts", "emb_test"); // TODO: replace absolute path
    mi.set_model_max_size(49917);
    int sz = mi.get_array_from_ctx_internal("vertex_positions", 0);
    mi.get_array_from_ctx_internal("vertex_normals", 1);
    mi.get_array_from_ctx_internal("vertex_texcoords", 2);
    int vertex_count = sz/3;
    logerr("wrote %d vertices (%d pos buffer size)", vertex_count, sz);
    std::vector<float> model_332 = std::vector<float>(8*vertex_count);
    for (int i=0;i<vertex_count;i++)
    {
      model_332[8*i] = mi.buffers[0][3*i];
      model_332[8*i+1] = mi.buffers[0][3*i+1];
      model_332[8*i+2] = mi.buffers[0][3*i+2];

      model_332[8*i+3] = mi.buffers[1][3*i];
      model_332[8*i+4] = mi.buffers[1][3*i+1];
      model_332[8*i+5] = mi.buffers[1][3*i+2];

      model_332[8*i+6] = mi.buffers[2][2*i];
      model_332[8*i+7] = mi.buffers[2][2*i+1];
    }
    for (int i = 0; i < 24; i++)
    {
      //logerr("%d model %f", i, model_332[i]);
    }

    size_t x_n = 4;
    std::vector<dgen::dfloat> X(x_n);
    std::vector<int> inds;
    std::vector<dgen::dfloat> Y;

    CppAD::Independent(X);
    {
      dgen::dmat43 sc = dgen::translate(dgen::rotate(dgen::ident(), dgen::dvec3{0,1,0}, X[0]), dgen::dvec3{X[1], X[2], X[3]});
      Y.reserve(model_332.size());
      for (float &f : model_332)
        Y.push_back(f);
      dgen::transform(Y, sc);
    }
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording

    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    
    X0[0] = 0;
    X0[1] = 0;
    X0[2] = 0.25;
    X0[3] = 0;

    for (int iter = 0; iter < 10; iter++)
    {
      jac = f.Jacobian(X0); // Jacobian for operation sequence
      //dgen::print_jackobian(jac, x_n, y_n, 100);
      res = f.Forward(0, X0); 
      for (int i = 0; i < 24; i++)
      {
        //logerr("%d res %f", i, res[i]);
      }

      for (int i=0;i<vertex_count;i++)
      {
        mi.buffers[0][3*i] = res[8*i];
        mi.buffers[0][3*i+1] = res[8*i+1];
        mi.buffers[0][3*i+2] = res[8*i+2];

        mi.buffers[1][3*i] = res[8*i+3];
        mi.buffers[1][3*i+1] = res[8*i+4];
        mi.buffers[1][3*i+2] = res[8*i+5];

        mi.buffers[2][2*i] = res[8*i+6];
        mi.buffers[2][2*i+1] = res[8*i+7];
      }

      mi.set_array_to_ctx_internal("vertex_positions", 0, 3*vertex_count);
      mi.set_array_to_ctx_internal("vertex_normals", 1, 3*vertex_count);
      mi.set_array_to_ctx_internal("vertex_texcoords", 2, 2*vertex_count);
      float loss = mi.render_and_compare_internal(4);

      mi.get_array_from_ctx_internal("vertex_positions_grad", 0);
      mi.get_array_from_ctx_internal("vertex_normals_grad", 1);
      mi.get_array_from_ctx_internal("vertex_texcoords_grad", 2);

      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      for (int i=0;i<vertex_count;i++)
      {
        for (int j=0;j<x_n;j++)
        {
          final_grad[j] += jac[(8*i)*x_n + j]*mi.buffers[0][3*i];
          final_grad[j] += jac[(8*i+1)*x_n + j]*mi.buffers[0][3*i+1];
          final_grad[j] += jac[(8*i+2)*x_n + j]*mi.buffers[0][3*i+2];

          final_grad[j] += jac[(8*i+3)*x_n + j]*mi.buffers[1][3*i];
          final_grad[j] += jac[(8*i+4)*x_n + j]*mi.buffers[1][3*i+1];
          final_grad[j] += jac[(8*i+5)*x_n + j]*mi.buffers[1][3*i+2];

          final_grad[j] += jac[(8*i+6)*x_n + j]*mi.buffers[2][2*i];
          final_grad[j] += jac[(8*i+7)*x_n + j]*mi.buffers[2][2*i+1];
        }
      }

      debug("[");
      for (int j=0;j<x_n;j++)
      {
        debug("%.6f ", X0[j]);
      }
      debug("]\n");
      debug("{");
      for (int j=0;j<x_n;j++)
      {
        debug("%.6f ", final_grad[j]);
        X0[j] = X0[j] - final_grad[j];
      }
      debug("]}\n");
    }
    mi.finish();
  }
}