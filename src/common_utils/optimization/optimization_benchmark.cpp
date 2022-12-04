#include "optimization.h"
#include "optimization_benchmark.h"
#include "common_utils/blk.h"

namespace opt
{
  void optimization_benchmark()
  {
      opt_func_with_grad F_sphere = [&](std::vector<float> &params) -> std::pair<float,std::vector<float>>
      {
        std::pair<float,std::vector<float>> res;
        res.first = 0;
        for (int i=0;i<params.size();i++)
        {
          res.first += (params[i]-i)*(params[i]-i);
          res.second.push_back(2*(params[i]-i));
        }
        return res;
      };

      opt_func_with_grad F_Rosenbrock = [&](std::vector<float> &params) -> std::pair<float,std::vector<float>>
      {
        std::pair<float,std::vector<float>> res;
        res.first = 0;
        for (int i=0;i<params.size() - 1;i++)
        {
          res.first += 100*SQR(params[i+1] - params[i]*params[i]) + SQR(1 - params[i]);
          if (i == 0)
            res.second.push_back(-400*params[0]*(params[1] - params[0]*params[0]) - 2*(1-params[0]));
          else
            res.second.push_back(200*(params[i] - params[i-1]*params[i-1])-400*params[i]*(params[i+1] - params[i]*params[i]) - 2*(1-params[i]));
        }
        res.second.push_back(200*(params[params.size() - 1] - params[params.size() - 2]*params[params.size() - 2]));
        return res;
      };

      opt_func_with_grad F_rastrigin = [&](std::vector<float> &params) -> std::pair<float,std::vector<float>>
      {
        std::pair<float,std::vector<float>> res;
        int N = params.size();
        int A = 10;
        res.first = A*N;
        for (int i=0;i<N;i++)
        {
          res.first += params[i]*params[i] - A*cos(2*PI*params[i]);
          res.second.push_back(2*params[i] + A*sin(2*PI*params[i]));
        }
        return res;
      };
    Optimizer *gd = new Adam();
    std::vector<float> min_v = {-5, -5, -5, -5, -5};
    std::vector<float> max_v = {4, 4, 4, 4, 4};
    Block settings;
    settings.set_int("iterations", 2500);
    settings.set_double("learning_rate", 20);
    gd->optimize(F_rastrigin, min_v, max_v, settings);
    gd->print_stat();
  }
}