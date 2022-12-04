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
  class Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings) = 0;
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
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings) override;
  };
  class Adam : public Optimizer
  {
  public:
    virtual void optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings) override;
  private:
    std::vector<float> V; 
    std::vector<float> S; 
    float alpha = 1;
    float beta_1 = 1;
    float beta_2 = 1;
    float eps = 1;
  };
}