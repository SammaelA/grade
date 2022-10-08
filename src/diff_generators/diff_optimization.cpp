#include "diff_optimization.h"
#include "diff_geometry_generation.h"
#include "mitsuba_python_interaction.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"
#include <functional>
#include <chrono>
#include "common_utils/blk.h"
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
    void init(dgen::generator_func _model_creator, int _model_creator_params_cnt)
    {
      model_creator = _model_creator;
      model_creator_params_cnt = _model_creator_params_cnt;
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
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        model_creator(X, Y);
        dgen::transform_by_scene_parameters(X, model_creator_params_cnt, Y);
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
    int model_creator_params_cnt = 0;
    dgen::generator_func model_creator;
    std::map<int, int> hash_to_function_pos;
    std::vector<CppAD::ADFun<float> *> functions;
    std::vector<int> output_sizes;//same size as functions vector
  };

  struct OptimizationUnitGD
  {
    void init(int _id, const std::vector<float> &init_params, DiffFunctionEvaluator &_func, MitsubaInterface &_mi, 
              const std::vector<float> &_params_min, const std::vector<float> &_params_max,
              bool _verbose = false)
    {
      assert(init_params.size() > 0);
      assert(init_params.size() == _params_min.size());
      assert(init_params.size() == _params_max.size());
      params_min = _params_min;
      params_max = _params_max;
      id = _id;
      verbose = _verbose;
      func = &_func;
      mi = &_mi;
      params = init_params;
      opt = new Adam2(0.015);
      x_n = init_params.size();
      for (int i=0;i<8;i++)
        timers[i] = 0;
    }
    void iterate()
    {
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    //float ms = 1e-4 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      std::vector<float> jac = func->get_jac(params);
      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      std::vector<float> res = func->get(params); 
      std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      float loss = mi->render_and_compare(res);
      std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
      mi->compute_final_grad(jac, x_n, res.size()/FLOAT_PER_VERTEX, final_grad);
      std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
      if (verbose)
      {
        debug("[%d] loss = %.3f, quality = %.3f\n", iterations, loss, get_quality());

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
      else if (iterations % 5 == 0)
        debug("[%d] loss = %.3f, quality = %.3f\n", iterations, loss, get_quality());
      std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();

      if (loss < best_error)
      {
        best_error = loss;
        best_params = params;
        best_error_iter = iterations;
      }
      best_error_stat.push_back(best_error);
      params = opt->step(params, final_grad, loss);
      for (int i=0;i<params.size();i++)
      {
        if (params[i] <= params_min[i])
        {
          float rnd = urand(0, 1);
          params[i] = params_min[i] + rnd*MIN(params_min[i] - params[i], params_max[i] - params_min[i]);
        }
        else if (params[i] >= params_max[i])
        {
          float rnd = urand(0, 1);
          params[i] = params_max[i] - rnd*MIN(params[i] - params_max[i], params_max[i] - params_min[i]);
        }
      }
      iterations++;

      std::chrono::steady_clock::time_point t8 = std::chrono::steady_clock::now();
      timers[0] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      timers[1] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
      timers[2] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t5 - t3).count();
      timers[3] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
      timers[4] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
      timers[5] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
    }

    ~OptimizationUnitGD()
    {
      if (opt)
        delete opt;
    }

    float get_quality()
    {
      float time_stat = (iterations <= 10) ? (11 - iterations) : (1 - 0.01*iterations);
      float val_stat = 1/(0.001 + best_error);
      float change_stat = 0;
      if (iterations >= 10)
        change_stat = 1 - best_error_stat[best_error_stat.size()-1 - 10] / (1e-9 + best_error);
      float decay_stat = MAX(0, 0.1*(iterations - best_error_iter - 10));
      return val_stat*exp(time_stat)*exp(change_stat)*exp(-decay_stat);
    }

    void print_stat()
    {
      float total_time = 1e-6;
      for (int i=0;i<8;i++)
        total_time += timers[i];
      total_time *= 1e-3;

      debug("Optimization Unit %d statistics\n", id);
      debug("Input parameters count: %d\n", x_n);
      debug("Itarations: %d\n", iterations);
      debug("Best value: %.4f\n", best_error);
      debug("Total time: %.3f s\n", total_time);
      std::vector<std::string> markers = {
        "Model jacobian calculation", "Model calculation", "Rendering", "Final Gradient calculation",
        "Debug print", "Optimizaiton"
      };
      for (int i=0;i<markers.size();i++)
      {
        float time = 1e-3*timers[i];
        float part = 100*time/total_time;
        float per_iter = 1000*time/iterations;
        debug("%s: %.3f s (%.1f %), %.1f ms per iteration\n", markers[i].c_str(), time, part, per_iter);
      }
    }

    DiffFunctionEvaluator *func = nullptr;
    MitsubaInterface *mi = nullptr;
    Optimizer *opt = nullptr;
    std::vector<float> params;
    std::vector<float> best_params;
    std::vector<float> params_min;
    std::vector<float> params_max;
    float best_error = 1e9;
    std::vector<float> best_error_stat;
    int best_error_iter = 0;
    int iterations = 0;
    int x_n = 0;
    bool verbose = false;
    double timers[8];
    int id = 0;
  };

  void test()
  {
    Block gen_params, scene_params;
    std::vector<float> params_min, params_max;
    load_block_from_file("dishes_gen_parameters_description.blk", gen_params);
    load_block_from_file("diff_gen_scene_parameters_description.blk", scene_params);

    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();
    auto process_blk = [&](Block &blk){
      for (int i=0;i<blk.size();i++)
      {
        Block *pb = blk.get_block(i);
        if (!pb && pb->size()>0)
        {
          logerr("invalid parameter description\"%s\". It should be non-empty block", blk.get_name(i));
        }
        else
        {
          glm::vec2 min_max = pb->get_vec2("values", glm::vec2(1e9,-1e9));
          if (min_max.x > min_max.y)
            logerr("invalid parameter description\"%s\". It should have values:p2 with min and max values", blk.get_name(i));
          params_min.push_back(min_max.x);
          params_max.push_back(min_max.y);
        }
      }
    };
    process_blk(gen_params);
    process_blk(scene_params);

    size_t x_n = gen_params_cnt + scene_params_cnt;
    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene)\n", 
          x_n, gen_params_cnt, scene_params_cnt);
    
    std::vector<float> reference_params{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                                        0.08, 0.25, 0.5, //hand params
                                        0, PI/4, 0, 0, 0, 0};//rotation and transform
    std::vector<float> init_params{4, 4, 4, 4, 4, 4, 4, 4, 4,
                                   0.05, 0.3, 0.4,
                                   0, PI/5, 0, 0, 0, 0};
    DiffFunctionEvaluator func;
    func.init(dgen::create_cup, gen_params_cnt);

    std::vector<float> reference = func.get(reference_params);

    MitsubaInterface mi;
    mi.init("scripts", "emb_test", MitsubaInterface::RenderSettings(196, 196, 1, MitsubaInterface::MitsubaVariant::LLVM));
    mi.init_optimization("saves/reference.png", MitsubaInterface::LOSS_MSE, 1 << 16, false);
    mi.render_model_to_file(reference, "saves/reference.png");

    OptimizationUnitGD opt_unit;
    opt_unit.init(0, init_params, func, mi, params_min, params_max, false);
    for (int j=0;j<250;j++)
        opt_unit.iterate();
    opt_unit.print_stat();
    mi.finish();
  }
}