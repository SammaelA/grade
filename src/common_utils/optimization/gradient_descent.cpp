#include "optimization.h"
#include "common_utils/blk.h"
namespace opt
{
  void GradientDescent::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                                 init_params_func &get_init_params)
  {
    assert(min_X.size() > 0);
    assert(min_X.size() == max_X.size());

    int iterations = settings.get_int("iterations");
    float lr = settings.get_double("learning_rate");
    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;

    int N = min_X.size();
    std::vector<float> X;
    settings.get_arr("initial_params", X);
    if (X.empty())
    {
      for (int i=0;i<N;i++)
        X.push_back(0.5*(min_X[i] + max_X[i]));
    }
    assert(X.size() == min_X.size());
    
    for (int iter=0;iter<iterations;iter++)
    {
      std::vector<std::vector<float>> Xv = {X};
      auto res = F(Xv);
      float &val = res[0].first;
      std::vector<float> &grad = res[0].second;
      if (val < best_result)
      {
        best_params = X;
        best_result = val;
      }
      for (int i=0;i<X.size();i++)
        X[i] = CLAMP(X[i] - lr*grad[i], min_X[i], max_X[i]);
      if ((iter % 100 == 0) && verbose)
        debug("Adam iter %d val = %.4f best_val = %.4f\n", iter, val, best_result);
    }
  }

  void Adam::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                      init_params_func &get_init_params)
  {
    assert(min_X.size() > 0);
    assert(min_X.size() == max_X.size());

    int iterations = settings.get_int("iterations");

    float alpha = settings.get_double("learning_rate", 10);
    float beta_1 = settings.get_double("beta_1", 0.9);
    float beta_2 = settings.get_double("beta_2", 0.999);
    float eps = settings.get_double("eps", 1e-8);
    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;

    int N = min_X.size();
    std::vector<float> X;
    settings.get_arr("initial_params", X);
    if (X.empty())
    {
      for (int i=0;i<N;i++)
        X.push_back(0.5*(min_X[i] + max_X[i]));
    }
    assert(X.size() == min_X.size());

    std::vector<float> grad_mult;//this additional array can be used as mask to prevent changing some of parameters 
    settings.get_arr("derivatives_mult", grad_mult);
    if (grad_mult.empty())
      grad_mult = std::vector<float>(N, 1);//default value - change all parameters
    assert(grad_mult.size() == min_X.size());

    std::vector<float> V = std::vector<float>(X.size(), 0); 
    std::vector<float> S = std::vector<float>(X.size(), 0);

    for (int iter=0; iter<iterations; iter++)
    {
      std::vector<std::vector<float>> Xv = {X};
      auto res = F(Xv);
      float &val = res[0].first;
      std::vector<float> &x_grad = res[0].second;
      if (val < best_result)
      {
        best_params = X;
        best_result = val;
      }
      for (int i=0;i<X.size();i++)
      {
        float g = grad_mult[i]*x_grad[i];
        V[i] = beta_1 * V[i] + (1-beta_1)*g;
        float Vh = V[i] / (1 - pow(beta_1, iter+1)); 
        S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
        float Sh = S[i] / (1 - pow(beta_2, iter+1)); 
        X[i] = CLAMP(X[i] - alpha*Vh/(sqrt(Sh) + eps), min_X[i], max_X[i]);
      }
      if ((iter % 5 == 0) && verbose)
        debug("Adam iter %d val = %.4f best_val = %.4f\n", iter, val, best_result);
    }
  }
}