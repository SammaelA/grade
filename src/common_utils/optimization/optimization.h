#pragma once
#include <vector>
#include <string>
#include <functional>
#include "common_utils/utility.h"
struct Block;

namespace opt
{
  typedef std::function<float(std::vector<float> &)> opt_func;
  typedef std::function<std::pair<float,std::vector<float>>(std::vector<float> &)> opt_func_with_grad;
  typedef std::function<std::vector<float>(std::vector<std::vector<float>> &)> opt_func_vector;
  typedef std::function<std::vector<std::pair<float,std::vector<float>>>(std::vector<std::vector<float>> &)> opt_func_with_grad_vector;

  void check_gradients(opt_func_with_grad &F, const std::vector<float> &min_X, const std::vector<float> &max_X,
                       int samples = 100, float h = 0.001, bool show_detailed_info = false);

  class Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) = 0;
    virtual std::vector<float> get_best_result(float *val = nullptr)
    {
      if (val)
        *val = best_result;
      return best_params;
    }
    Optimizer() = default;
    virtual ~Optimizer() {};
    virtual void get_n_best(int n, std::vector<std::vector<float>> &params, std::vector<float> *values = nullptr)
    {
      logerr("This saves only one best result. Use get_best_result() to get it");
    }
    virtual void print_stat()
    {
      float val = -1;
      std::vector<float> params = get_best_result(&val); 
      debug("Optimization finished\n");
      debug("Best parameters: (");
      for (float &v : params)
      {
        debug("%.4f, ", v);
      }
      debug(")\n");
      debug("Best value: %.5f\n",val);
    }

    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings)
    {
      opt_func_vector F_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<float>
      {
        std::vector<float> res = std::vector<float>(params.size(), 0);
        return res;
      };
      optimize(F, min_X, max_X, settings, F_vec);
    }

    virtual void optimize(opt_func_with_grad &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func &F_reg)
    {
      opt_func_with_grad_vector F_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<std::pair<float,std::vector<float>>>
      {
        std::vector<std::pair<float,std::vector<float>>> res;
        for (auto &x : params)
          res.push_back(F(x));
        return res;
      };
      opt_func_vector F_reg_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<float>
      {
        std::vector<float> res;
        for (auto &x : params)
          res.push_back(F_reg(x));
        return res;
      };
      optimize(F_vec, min_X, max_X, settings, F_reg_vec);
    }

    virtual void optimize(opt_func &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings)
    {
      opt_func_with_grad_vector F_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<std::pair<float,std::vector<float>>>
      {
        std::vector<std::pair<float,std::vector<float>>> res;
        for (auto &x : params)
          res.push_back(std::pair<float,std::vector<float>>(F(x), std::vector<float>()));
        return res;
      };
      optimize(F_vec, min_X, max_X, settings);
    }

    virtual void optimize(opt_func_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings)
    {
      opt_func_with_grad_vector F_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<std::pair<float,std::vector<float>>>
      {
        std::vector<std::pair<float,std::vector<float>>> res;
        std::vector<float> vals = F(params);
        for (auto &v : vals)
          res.push_back(std::pair<float,std::vector<float>>(v, std::vector<float>()));
        return res;
      };
      optimize(F_vec, min_X, max_X, settings);
    }

    virtual void optimize(opt_func_with_grad &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings)
    {
      opt_func_with_grad_vector F_vec = [&](std::vector<std::vector<float>> &params) -> std::vector<std::pair<float,std::vector<float>>>
      {
        std::vector<std::pair<float,std::vector<float>>> res;
        for (auto &x : params)
          res.push_back(F(x));
        return res;
      };
      optimize(F_vec, min_X, max_X, settings);
    }
  protected:
    std::vector<float> best_params;
    float best_result = 1e9;
  };

  class GradientDescent : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  };
  class Adam : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  };
  class DifferentialEvolutionOptimizer : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  } ;

  class GridSearchAdam : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  } ;
  class CMA_ES : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  } ;
  class MemeticClassic : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                          opt_func_vector &F_reg) override;
  } ;
}