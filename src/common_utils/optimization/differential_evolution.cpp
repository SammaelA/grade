#include "optimization.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"

namespace opt
{
  void DifferentialEvolutionOptimizer::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, 
                                                Block &settings, opt_func_vector &F_reg)
  {
    /*
    Storn R., Price K. Differential evolution–a simple and efficient heuristic for global optimization over continuous spaces
    //Journal of global optimization. – 1997. – Т. 11. – №. 4. – С. 341-359.
    */

    assert(min_X.size() > 0);
    assert(min_X.size() == max_X.size());

    int generations = settings.get_int("generations", 10);
    int population_size = settings.get_int("population_size", 10);
    float K = settings.get_double("K", 0.7);
    float crossover_rate = settings.get_double("crossover_rate", 0.3);
    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;
    int local_search_iterations = settings.get_int("local_search_iterations", 40);
    float local_search_learning_rate = settings.get_double("local_search_learning_rate", 0.5);

    int N = min_X.size();

    std::vector<std::vector<float>> population(population_size, std::vector<float>(N,0));
    std::vector<float> values(population_size, 1e9);

    std::vector<std::vector<float>> new_population(population_size, std::vector<float>(N,0));
    std::vector<float> new_values(population_size, 1e9);

    for (int i=0;i<population_size;i++)
    {
      for (int j=0;j<N;j++)
        population[i][j] = urand(min_X[j], max_X[j]);
    }
    {
      auto results = F(population);
      for (int i=0;i<population_size;i++)
        values[i] = results[i].first;
    }
    for (int generation = 0; generation < generations; generation++)
    {
      for (int i=0;i<population_size;i++)
      {
        auto X_r = population[rand() % population_size];
        auto X_s = population[rand() % population_size];
        auto X_t = population[rand() % population_size];

        for (int j=0;j<N;j++)
        {
          if (urand() < crossover_rate)
            new_population[i][j] = population[i][j];
          else
            new_population[i][j] = X_t[j] + urand(0,K)*(X_r[j] - X_s[j]);
        }
        Optimizer *opt = new Adam();
        Block adam_settings;
        adam_settings.add_arr("initial_params", new_population[i]);
        adam_settings.add_double("learning_rate", local_search_learning_rate);
        adam_settings.add_int("iterations", local_search_iterations);
        opt->optimize(F, min_X, max_X, adam_settings);
        new_population[i] = opt->get_best_result();
        delete opt;
      }

      auto new_results = F(new_population);

      for (int i=0;i<population_size;i++)
      {
        if (new_results[i].first < values[i])
        {
          population[i] = new_population[i];
          values[i] = new_results[i].first;
        }
        if (values[i] < best_result)
        {
          best_result = values[i];
          best_params = population[i];
        }
      }
      if (verbose)
        debug("DE generation %d. Best value %.4f\n", generation, best_result);
    }
  }
}