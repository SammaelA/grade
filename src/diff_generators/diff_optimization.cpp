#include "diff_optimization.h"
#include "diff_geometry_generation.h"
#include "mitsuba_python_interaction.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"
#include <functional>

namespace dopt
{
  class Optimizer
  {
  public:
    Optimizer(){};
    virtual ~Optimizer(){};
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) = 0;
  };

  class GradientDescentSimple : public Optimizer
  {
  public:
    GradientDescentSimple(float _alpha = 0.1) 
    {
      assert(_alpha > 0);
      alpha = _alpha;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
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
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
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
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
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

  class Adam2 : public Optimizer
  {
  public:
    Adam2(float _alpha = 0.01, float _beta_1 = 0.9, float _beta_2 = 0.999, float _eps = 1e-8)
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
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
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
      std::vector<float> x(x_prev.size(),0);

      if (value < 1.2*prev_val)
      {
        iter++;
        for (int i=0;i<x_prev.size();i++)
        {
          V[i] = beta_1 * V[i] + (1-beta_1)*x_grad[i];
          float Vh = V[i] / (1 - pow(beta_1, iter)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*x_grad[i]*x_grad[i];
          float Sh = S[i] / (1 - pow(beta_2, iter)); 

          x[i] = x_prev[i] - alpha*Vh/(sqrt(Sh) + eps);
        }
        prev_x = x_prev; 
        prev_val = value;
      }
      else
      {
        for (int i=0;i<x_prev.size();i++)
        {
          float rnd = urand();
          x[i] = rnd*prev_x[i] + (1-rnd)*x_prev[i];
        } 
      }
      return x;
    }
  private:
    std::vector<float> V; 
    std::vector<float> S; 
    std::vector<float> prev_x;
    float prev_val = 1e9;
    float alpha = 1;
    float beta_1 = 1;
    float beta_2 = 1;
    float eps = 1;
    int iter = 0;
  };

  class DiffFunctionEvaluator
  {
  public:
    ~DiffFunctionEvaluator()
    {
      for (auto f : functions)
      {
        if (f)
          delete f;
      }
    }
    void init(std::function<void(std::vector<dgen::dfloat> &X, std::vector<dgen::dfloat> &Y)> _model_creator)
    {
      model_creator = _model_creator;
    }
    std::vector<float> get(const std::vector<float> &params)
    {
      return functions[find_or_add(params)]->Forward(0, params); 
    }
    std::vector<float> get_jac(const std::vector<float> &params)
    {
      return functions[find_or_add(params)]->Jacobian(params); 
    }
  private:
    int find_or_add(const std::vector<float> &params)
    {
      int hash = get_function_hash(params);
      auto it = hash_to_function_pos.find(hash);
      if (it == hash_to_function_pos.end())
      {
        std::vector<dgen::dfloat> X(params.size());
        //std::vector<float> test_params = {3.264063, 3.752163, 4.444632, 4.231698, 3.847324, 3.966257, 3.936437, 4.553413, 4.441343, 0.050000, 0.441019, 0.468129, 3.219785};
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        model_creator(X, Y);
        CppAD::ADFun<float> *f = new CppAD::ADFun<float>(X, Y); 
        functions.push_back(f);
        output_sizes.push_back(Y.size());
        int f_pos = functions.size()-1;
        hash_to_function_pos.emplace(hash, f_pos);

        return f_pos;
      }
      return it->second;
    }

    int get_function_hash(const std::vector<float> &params)
    {
      return 0;
    }
    std::function<void(std::vector<dgen::dfloat> &X, std::vector<dgen::dfloat> &Y)> model_creator;
    std::map<int, int> hash_to_function_pos;
    std::vector<CppAD::ADFun<float> *> functions;
    std::vector<int> output_sizes;//same size as functions vector
  };

  struct OptimizationUnitGD
  {
    void init(const std::vector<float> &init_params, DiffFunctionEvaluator &_func, MitsubaInterface &_mi, bool _verbose = false)
    {
      verbose = _verbose;
      func = &_func;
      mi = &_mi;
      params = init_params;
      opt = new Adam2(0.015);
      x_n = init_params.size();
    }
    void iterate()
    {
      std::vector<float> jac = func->get_jac(params);
      std::vector<float> res = func->get(params); 
      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      float loss = mi->render_and_compare(res);
      mi->compute_final_grad(jac, x_n, res.size()/FLOAT_PER_VERTEX, final_grad);

      if (verbose)
      {
        debug("[%d] loss = %.3f\n", iterations, loss);

        debug("params [");
        for (int j=0;j<x_n;j++)
        {
          debug("%.3f, ", params[j]);
        }
        debug("]\n");

        debug("grad {");
        for (int j=0;j<x_n;j++)
        {
          debug("%.3f ", final_grad[j]);
        }
        debug("}\n");
      }
      else if (iterations % 10 == 0)
        debug("[%d] loss = %.3f\n", iterations, loss);

      iterations++;
      if (loss < best_error)
      {
        best_error = loss;
        best_params = params;
      }
      params = opt->step(params, final_grad, loss);
      for (int i=0;i<12;i++)
      {
        float &p = params[i];
        if (p > 5)
          p = 5;
        if (p < 0.05)
          p = 0.05;
      }
      float d = params[11] - params[10];
      if (d < 3 * params[9])
      {
        params[10] -= 1.5*params[9];
        params[11] += 1.5*params[9];
      }
    }
    ~OptimizationUnitGD()
    {
      if (opt)
        delete opt;
    }

    DiffFunctionEvaluator *func = nullptr;
    MitsubaInterface *mi = nullptr;
    Optimizer *opt = nullptr;
    std::vector<float> params;
    std::vector<float> best_params;
    float best_error = 1e9;
    int iterations = 0;
    int x_n = 0;
    bool verbose = false;
  };

  void test()
  {
    constexpr size_t x_n = 18;
    std::vector<float> reference_params{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                                        0.08, 0.17, 0.83, //hand params
                                        0, PI/4, 0, 0, 0, 0};//rotation and transform
    std::vector<float> init_params{4, 4, 4, 4, 4, 4, 4, 4, 4,
                                   0.05, 0.3, 0.7,
                                   0, PI/5, 0, 0, 0, 0};
    
    DiffFunctionEvaluator func;
    func.init(dgen::create_cup);

    std::vector<float> reference = func.get(reference_params);

    MitsubaInterface mi;
    mi.init("scripts", "emb_test", MitsubaInterface::RenderSettings(196, 196, 1, MitsubaInterface::MitsubaVariant::LLVM));
    mi.init_optimization("saves/reference.png", MitsubaInterface::LOSS_MSE, 1 << 16);
    mi.render_model_to_file(reference, "saves/reference.png");

    OptimizationUnitGD opt_unit;
    opt_unit.init(init_params, func, mi, false);
    for (int j=0;j<20;j++)
        opt_unit.iterate();

    mi.finish();
  }
}