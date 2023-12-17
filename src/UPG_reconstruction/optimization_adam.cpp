#include "optimization.h"

namespace upg
{
  class UPGOptimizerAdam : public UPGOptimizer
  {
  public:
    UPGOptimizerAdam(UPGOptimizableFunction *_func, 
                     const Block &settings, const UPGReconstructionResult &start_params) :
    UPGOptimizer(_func)
    {
      gen = func->get_generator(start_params.structure);
      iterations = settings.get_int("iterations");
      alpha = settings.get_double("learning_rate", 0.01);
      beta_1 = settings.get_double("beta_1", 0.9);
      beta_2 = settings.get_double("beta_2", 0.999);
      eps = settings.get_double("eps", 1e-8);
      verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;

      pd = func->get_full_parameters_description(gen.get());
      X_n = start_params.parameters.p.size();
      X = func->gen_params_to_opt_params(start_params.parameters.p, pd);
      gen_structure = start_params.structure;


      V = std::vector<float>(X_n, 0); 
      S = std::vector<float>(X_n, 0);
      best_params = X;
      x_grad = std::vector<float>(X_n, 0); 
      best_result = 1e9;
    }
    virtual void optimize(int iters = -1) override
    {
      for (int iter=0; iter< (iters>0 ? iters : iterations); iter++)
      {
        float val = func->f_grad_f(gen.get(), pd, X, x_grad);
        if (val < best_result)
        {
          best_params = X;
          best_result = val;
        }
        for (int i=0;i<X_n;i++)
        {
          float g = x_grad[i];
          V[i] = beta_1 * V[i] + (1-beta_1)*g;
          float Vh = V[i] / (1 - pow(beta_1, total_iterations+1)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
          float Sh = S[i] / (1 - pow(beta_2, total_iterations+1)); 
          X[i] -= alpha*Vh/(sqrt(Sh) + eps);
        }
        if ((total_iterations % 5 == 0) && verbose)
          debug("Adam iter %3d  val = %.8f best_val = %.8f\n", total_iterations, val, best_result);
        total_iterations++;
      }
      if (verbose)
        debug("Adam final res val = %.8f best_val = %.8f\n", best_result, best_result);
    }

    virtual std::vector<UPGReconstructionResult> get_best_results() override
    {
      UPGReconstructionResult res;
      res.structure = gen_structure;
      {
        res.parameters.p = func->opt_params_to_gen_params(best_params, pd);
        res.loss_optimizer = best_result;
      }
      return {res};
    }
    OptParams get_best_params_in_optimizer_format() { return best_params;}
    float get_best_result_in_optimizer_format() { return best_result;}

  private:
    int iterations;
    float alpha, beta_1, beta_2, eps;
    bool verbose;
    int X_n;
    OptParams X;
    std::shared_ptr<UniversalGenInstance> gen;
    ParametersDescription pd;
    UPGStructure gen_structure;

    std::vector<float> V; 
    std::vector<float> S;
    OptParams best_params;
    std::vector<float> x_grad; 
    float best_result;
    int total_iterations = 0;
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_adam(UPGOptimizableFunction *_func, 
                                                   const Block &settings, const UPGReconstructionResult &start_params)
  {
    return std::make_shared<UPGOptimizerAdam>(_func, settings, start_params);
  }
}