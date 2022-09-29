#include "diff_optimization.h"
#include "diff_geometry_generation.h"
#include "mitsuba_python_interaction.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"

namespace dopt
{
  class Optimizer
  {
  public:
    Optimizer(){};
    virtual ~Optimizer(){};
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad) = 0;
  };

  class GradientDescentSimple : public Optimizer
  {
  public:
    GradientDescentSimple(float _alpha = 0.1) 
    {
      assert(_alpha > 0);
      alpha = _alpha;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad) override
    {
      assert(x_prev.size() == x_grad.size());
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        x[i] = x_prev[i] - alpha*x_grad[i];
      }

      return x;
    }
  private:
    float alpha = 1;
  };

  class RMSprop : public Optimizer
  {
  public:
    RMSprop(float _alpha = 0.01, float _beta = 0.99, float _eps = 1e-8)
    {
      assert(_alpha > 0);
      assert(_beta > 0);
      assert(_beta < 1);
      assert(_eps > 0);

      alpha = _alpha;
      beta = _beta;
      eps = _eps;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad) override
    {
      assert(x_prev.size() == x_grad.size());
      if (S.empty())
        S = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == S.size());
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        S[i] = beta * S[i] + (1-beta)*x_grad[i]*x_grad[i];
        x[i] = x_prev[i] - alpha*x_grad[i]/(sqrt(S[i]) + eps);
      }

      return x;
    }
  private:
    std::vector<float> S; 
    float alpha = 1;
    float beta = 1;
    float eps = 1;
  };

  class Adam : public Optimizer
  {
  public:
    Adam(float _alpha = 0.01, float _beta_1 = 0.9, float _beta_2 = 0.999, float _eps = 1e-8)
    {
      assert(_alpha > 0);
      assert(_beta_1 > 0);
      assert(_beta_1 < 1);
      assert(_beta_2 > 0);
      assert(_beta_2 < 1);
      assert(_eps > 0);
      
      alpha = _alpha;
      beta_1 = _beta_1;
      beta_2 = _beta_2;
      eps = _eps;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad) override
    {
      assert(x_prev.size() == x_grad.size());
      if (S.empty())
        S = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == S.size());
      if (V.empty())
        V = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == V.size());

      iter++;
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        V[i] = beta_1 * V[i] + (1-beta_1)*x_grad[i];
        float Vh = V[i] / (1 - pow(beta_1, iter)); 
        S[i] = beta_2 * S[i] + (1-beta_2)*x_grad[i]*x_grad[i];
        float Sh = S[i] / (1 - pow(beta_2, iter)); 

        x[i] = x_prev[i] - alpha*Vh/(sqrt(Sh) + eps);
      }

      return x;
    }
  private:
    std::vector<float> V; 
    std::vector<float> S; 
    float alpha = 1;
    float beta_1 = 1;
    float beta_2 = 1;
    float eps = 1;
    int iter = 0;
  };

  void test()
  {    
    size_t x_n = 3;
    std::vector<dgen::dfloat> X(x_n);
    int model_size = 0;
    std::vector<dgen::dfloat> Y;

    CppAD::Independent(X);
    {
      dgen::create_cup(X, Y);
    }
    size_t y_n = Y.size();
    int vertex_count = y_n/FLOAT_PER_VERTEX;
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording

    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    
    X0[0] = 4 - 1.45;
    //X0[1] = 4 - 1.0;
    //X0[2] = 4 - 0.65;
    //X0[3] = 4 - 0.45;
    X0[1] = 4 - 0.25;
    //X0[5] = 4 - 0.18;
    //X0[6] = 4 - 0.1;
    //X0[7] = 4 - 0.05;
    X0[2] = 4 - 0;

    res = f.Forward(0, X0); 
    
    MitsubaInterface mi;
    mi.init("scripts", "emb_test");
    mi.init_optimization("saves/reference.png", MitsubaInterface::RenderSettings(128, 128, 64), MitsubaInterface::LOSS_MSE_SQRT, 1 << 16);
    mi.render_model_to_file(res, MitsubaInterface::RenderSettings(512, 512, 64), "saves/reference.png");

    X0[0] = 4;
    X0[1] = 4;
    X0[2] = 4;
    //X0[3] = 4;
    //X0[4] = 4;
    //X0[5] = 4;
    //X0[6] = 4;
    //X0[7] = 4;
    //X0[8] = 4;

    Optimizer *opt = new Adam(0.025);

    int steps = 100;
    for (int iter = 0; iter < steps; iter++)
    {
      debug("[");
      for (int j=0;j<x_n;j++)
      {
        debug("%.6f ", X0[j]);
      }
      debug("]\n");

      jac = f.Jacobian(X0);
      //dgen::print_jackobian(jac, x_n, y_n, 100);
      res = f.Forward(0, X0); 
      for (int i = 0; i < 24; i++)
      {
        //logerr("%d res %f", i, res[i]);
      }
      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      float loss = mi.render_and_compare(res);
      debug("[%d/%d] loss = %.4f\n", iter, steps, loss);
      /*dgen::print_model(res);
      dgen::print_jackobian(jac, vertex_count, vertex_count, 100000);
      for (int i = 0; i < vertex_count; i++)
      {
        debug("%f %f %f\n", mi.buffers[0][3 * i], mi.buffers[0][3 * i+1], mi.buffers[0][3 * i+2]);
      }*/
      mi.compute_final_grad(jac, x_n, vertex_count, final_grad);

      debug("{");
      for (int j=0;j<x_n;j++)
      {
        debug("%.6f ", final_grad[j]);
      }
      debug("}\n");
      X0 = opt->step(X0, final_grad);
    }

    delete opt;
    mi.finish();
  }
}