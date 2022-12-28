#include "optimization.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"
#include <algorithm>

namespace opt
{
  void MemeticClassic::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                                opt_func_vector &F_reg)
  {
    float mutation_chance = settings.get_double("mutation_chance", 0.2);
    float mutation_power = settings.get_double("mutation_power", 0.3);
    int tournament_size = settings.get_int("tournament_size", 16);
    int population_size = settings.get_int("population_size", 32);
    int total_function_calls = settings.get_int("total_function_calls", 4800);
    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;
    int local_search_iterations = settings.get_int("local_search_iterations", 76);
    float local_search_learning_rate = settings.get_double("local_search_learning_rate", 0.05);
    float recreation_diversity_thr = settings.get_double("recreation_diversity_thr", 1.75);
    float depth_reg_q = settings.get_double("depth_reg_q", 0);
    int budget = 0;

    auto mutate = [&](const std::vector<float> &base) -> std::vector<float>
    {
      std::vector<float> res = base;
      for (int i=0; i<base.size(); i++)
      {
        if (urand() < mutation_chance)
          res[i] = CLAMP(base[i] + urand(-mutation_power, mutation_power)*(max_X[i] - min_X[i]), min_X[i], max_X[i]);
      }
      return res;
    };

    auto one_dot_crossover = [&](const std::vector<float> &A, const std::vector<float> &B) -> std::vector<float>
    {
      std::vector<float> res = A;

      int pos_st = urandi(0, res.size());
      int pos_en = res.size();
      if (urand() > 0.5)
      {
        pos_en = pos_st;
        pos_st = 0; 
      }
      for (int i=pos_st; i<pos_en; i++)
        res[i] = B[i];
      
      return res;
    };

    auto choose_parents_tournament = [&](const std::vector<std::vector<float>> &population, const std::vector<float> &values) -> std::pair<int, int>
    {
      std::vector<int> indexes = std::vector<int>(tournament_size, 0);
      for (int i=0;i<tournament_size;i++)
        indexes[i] = urandi(0, population.size());
      std::sort(indexes.begin(), indexes.end(), [&](const int & a, const int & b) -> bool{    
            return values[a] < values[b];});

      return std::pair<int, int>(indexes[0], indexes[1]);
    };

    auto initialize_population = [&](std::vector<std::vector<float>> &population)
    {
      population = std::vector<std::vector<float>>(population_size, std::vector<float>(min_X.size(), 0));
      for (int i=0;i<population_size;i++)
      {
        for (int j=0;j<min_X.size();j++)
          population[i][j] = urand(min_X[j], max_X[j]);
      }
    };

    auto local_search = [&](std::vector<float> start_params) -> std::pair<float, std::vector<float>>
    {
      Optimizer *opt = new Adam();
      Block adam_settings;
      adam_settings.add_arr("initial_params", start_params);
      adam_settings.add_double("learning_rate", local_search_learning_rate);
      adam_settings.add_int("iterations", local_search_iterations);
      adam_settings.add_bool("verbose", false);
      opt->optimize(F, min_X, max_X, adam_settings);
      budget += local_search_iterations;

      std::pair<float, std::vector<float>> res;
      res.second = opt->get_best_result(&(res.first));
      //add regualized part to the last result
      if (depth_reg_q > 0)
      {
        std::vector<std::vector<float>> t{res.second};
        res.first += depth_reg_q*F_reg(t)[0];
      }
      if (res.first < best_result)
      {
        best_result = res.first;
        best_params = res.second;
        logerr("%d new best %.4f", budget, best_result);
      }
      delete opt;
      return res;
    };

    auto calc_diversity = [&](std::vector<std::vector<float>> &population) -> float
    {
      float diversity = 0;
      int cnt = 0;
      for (int i1=0;i1<population.size();i1++)
      {
        for (int i2=i1+1;i2<population.size();i2++)
        {
          float d = 0;
          for (int j=0;j<min_X.size();j++)
            d += SQR((population[i1][j] - population[i2][j])/(max_X[j] - min_X[j]));
          diversity += sqrtf(d) / min_X.size();
          cnt++;
        }
      }
      diversity /= cnt;
      return diversity;
    };

    std::vector<std::vector<float>> population;
    std::vector<float> values;
    std::vector<float> qa_values;
    std::vector<int> local_improvements;
    std::vector<int> indices;

    initialize_population(population);
    for (int i=0; i<population.size(); i++)
    {
      auto res = local_search(population[i]);

      population[i] = res.second;
      values.push_back(res.first);
      qa_values.push_back(1e-3 / (res.first * res.first));
      local_improvements.push_back(1);
      indices.push_back(i);
    }
    float initial_diversity = calc_diversity(population);

    while (budget < total_function_calls)
    {
      //crossover + mutation
      auto parents = choose_parents_tournament(population, values);
      std::vector<float> new_solution = one_dot_crossover(population[parents.first], population[parents.second]);
      new_solution = mutate(new_solution);
      auto res = local_search(new_solution);

      int worst_idx = -1;
      float worst_val = res.first;
      for (int i=0;i<population.size();i++)
      {
        if (values[i] > worst_val)
        {
          worst_val = values[i];
          worst_idx = i;

        }
      }
      if (worst_idx >= 0)
      {
        logerr("replace with better child %.4f --> %.4f", values[worst_idx], res.first);
        values[worst_idx] = res.first;
        qa_values[worst_idx] = 1e-3 / (res.first * res.first);
        local_improvements[worst_idx] = 1;
        population[worst_idx] = res.second;
      }
      else
      {
        logerr("Child is the worst %.4f", res.first);
      } 

      //further local improvement
      if (urand() < ((float)budget)/total_function_calls)
      {
        float qav_sum = 0;
        for (float &v : qa_values)
          qav_sum += v;
        float rnd = urand(qav_sum);
        int improve_idx = 0;
        for (int i=0;i<population.size();i++)
        {
          if (rnd < qa_values[i])
          {
            improve_idx = i;
            break;
          }
          else
          {
            rnd -= qa_values[i];
          }
        }
        auto res = local_search(population[improve_idx]);
        logerr("improving existing solution (%d times) %.4f --> %.4f", local_improvements[improve_idx], values[improve_idx], res.first);
        population[improve_idx] = res.second;
        values[improve_idx] = res.first;
        local_improvements[improve_idx]++;
        //if we improve solution, there is unlikely that it will be further improved
        qa_values[improve_idx] = 1e-3 / (res.first * res.first * local_improvements[improve_idx]);
      }

      //diversity calculation and population restart
      float diversity = calc_diversity(population);
      logerr("population diversity %.3f", diversity);
      if (initial_diversity/(diversity + 1e-6) > recreation_diversity_thr)
      {
        //recreate population and place best solution in it
        logerr("recreating population");
        initialize_population(population);
        population[0] = best_params;

        for (int i=0; i<population.size(); i++)
        {
          auto res = local_search(population[i]);

          population[i] = res.second;
          values.push_back(res.first);
          qa_values.push_back(1e-3 / (res.first * res.first));
          local_improvements.push_back(1);
          indices.push_back(i);
        }
        initial_diversity = calc_diversity(population);        
      }
      
    }
  }
}