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


    int last_mutation_id = -1;
    std::array<int, 8> mutation_stat = {0,0,0,0,0,0,0,0}; //[10 9 2 16 7 8 0 1]
    auto mutate = [&](const std::vector<float> &base) -> std::vector<float>
    {
      last_mutation_id = -1;
      std::vector<float> res = base;
      for (int i=0; i<base.size(); i++)
      {
        if (urand() < mutation_chance)
          res[i] = CLAMP(base[i] + urand(-mutation_power, mutation_power)*(max_X[i] - min_X[i]), min_X[i], max_X[i]);
      }
      return res;
    };

    auto special_mutation = [&](const std::vector<float> &base) -> std::vector<float>
    {
      std::vector<float> res = base;
      int mutation_type = urandi(0, 5);
      last_mutation_id = mutation_type;
      if (mutation_type == 0)
      {
        float mul = urand(0.75, 1.25);
        for (int i=0; i<9; i++)
          res[i] = mul * base[i];
      }
      else if (mutation_type == 1)
      {
        res[10] = urand(0, 1);
      }
      else if (mutation_type == 2)
      {
        int end = base.size();
        float rnd = urand(0,1);
        if (rnd < 0.25)
          res[end - 1] += urand(-1, 1);
        else if (rnd < 0.75)
          res[end - 2] += urand(-1, 1);
        else
          res[end - 3] += urand(-1, 1);
      }
      else if (mutation_type == 3)
      {
        res[base.size() - 2] += urand(-3, 3);
      }
      else
      {
        for (int i=0; i<base.size(); i++)
        {
          if (urand() < mutation_chance)
            res[i] = CLAMP(base[i] + urand(-mutation_power, mutation_power)*(max_X[i] - min_X[i]), min_X[i], max_X[i]);
        }
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

      std::vector<std::vector<float>> presets_1 = {
        {2.014, 3.424, 3.404, 3.291, 3.276, 3.284, 3.357, 3.383, 3.354, 0.781, 1, 0.05, 0.7, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.4},
        {1.879, 1.888, 2.867, 3.070, 3.333, 3.533, 3.746, 3.899, 3.92, 0.690, 0.261, 0.055, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275},
        {5.5, 5.8, 6.1, 6.4, 6.7, 7, 7.3, 7.6, 7.9, 0.1, 0, 0.05, 0.7, 0.15, 0.25, 0.4, 0.35, 0.3, 0.25, 0.2},
        {2.014, 3.424, 3.404, 3.291, 3.276, 3.284, 3.357, 3.383, 3.354, 0.781, 1, 0.05, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3}
      };

      std::vector<std::vector<float>> presets_2 = {
        {0.1, 0.2, -0.2, 0.0, 0.0, 0.1},
        {0.1, 0, 0.004, 0.3, 3.141, -0.250},
        {0.1, 0, 0.004, 0.0, 0.0, 0.1},
        {0.1, 0.2, -0.2, 0.05, 1, 0.0},
        {0.0, 0.3, -0.3, 0.05, 0.5, 0.2},
        {0.1, 0.2, -0.2, 0.1, 2, 0.3},
        {0.0, 0.1, -0.25, 0.05, 1.5, 0.2},
        {0.1, 0.2, -0.2, 0.2, 0.0, 0.03}
      };

      std::vector<std::vector<float>> all_presets;
      for (auto &p1 : presets_1)
      {
        for (auto &p2 : presets_2)
        {
          all_presets.emplace_back();
          for (auto &v1 : p1)
            all_presets.back().push_back(v1);
          for (auto &v2 : p2)
            all_presets.back().push_back(v2);
        }
      }

      for (int i=0;i<population_size - 4;i++)
      {
        population[i] = all_presets[(int)urandi(0, all_presets.size())];
      }
    };

    auto local_search = [&](std::vector<float> start_params) -> std::pair<float, std::vector<float>>
    {
      Optimizer *opt = new Adam();
      Block adam_settings;
      int iterations = urandi(5, 2*local_search_iterations);
      adam_settings.add_arr("initial_params", start_params);
      adam_settings.add_double("learning_rate", local_search_learning_rate);
      adam_settings.add_int("iterations", iterations);
      adam_settings.add_bool("verbose", false);
      opt->optimize(F, min_X, max_X, adam_settings);
      budget += iterations;

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

    Block backup;

    if (load_block_from_file("backup.blk", backup))
    {
      bool cont = false;
      backup.get_int("continue", cont);
      if (cont)
      {
        backup.get_int("budget", budget);
        int sz = 0;
        backup.get_int("pop_size", sz);
        population.resize(sz);
        std::string num_pop = "AA_pop";
        for (int i = 0; i < population.size(); ++i)
        {
          num_pop[0] = 'A' + i / 26;
          num_pop[1] = 'A' + i % 26;
          backup.get_arr(num_pop.c_str(), population[i]);
        }
        backup.get_arr("values", values);
        backup.get_arr("qa_values", qa_values);
        backup.get_arr("local_improvements", local_improvements);
        backup.get_arr("indices", indices);
      }
    }

    while (budget < total_function_calls)
    {
      backup.clear();
      backup.add_bool("continue", true);
      backup.add_int("budget", budget);
      backup.add_int("pop_size", population.size());
      std::string num_pop = "AA_pop";
      for (int i = 0; i < population.size(); ++i)
      {
        num_pop[0] = 'A' + i / 26;
        num_pop[1] = 'A' + i % 26;
        backup.add_arr(num_pop.c_str(), population[i]);
      }
      backup.add_arr("values", values);
      backup.add_arr("qa_values", qa_values);
      backup.add_arr("local_improvements", local_improvements);
      backup.add_arr("indices", indices);
      save_block_to_file("backup.blk", backup);

      std::pair<float, std::vector<float>> res;
      if (urand() < 0.1)
      {
        //crossover + mutation
        auto parents = choose_parents_tournament(population, values);
        std::vector<float> new_solution = one_dot_crossover(population[parents.first], population[parents.second]);
        new_solution = special_mutation(new_solution);
        res = local_search(new_solution);
      }
      else
      {
        int id = urandi(0, population.size());
        std::vector<std::vector<float>> new_solution = {population[id]};
        new_solution[0] = special_mutation(new_solution[0]);
        auto fres = F(new_solution);
        res.first = fres[0].first;
        res.second = new_solution[0];
        logerr("mutated %f --> %f", values[id], res.first); 
      }
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

        if (last_mutation_id >= 0)
          mutation_stat[last_mutation_id]++;
        logerr("mutation stat [%d %d %d %d %d %d %d %d]", mutation_stat[0], mutation_stat[1], mutation_stat[2], mutation_stat[3],
                                                          mutation_stat[4], mutation_stat[5], mutation_stat[6], mutation_stat[7]);
      }
      else
      {
        logerr("Child is the worst %.4f", res.first);
      } 

      //further local improvement
      if (urand() < ((float)budget)/total_function_calls)
      {
        float qav_sum = 0;
        debug("qa_values [");
        for (float &v : qa_values)
        {
          debug("%f ", v);
          qav_sum += v;
        }
        debugnl();
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
        values.clear();
        qa_values.clear();
        local_improvements.clear();
        indices.clear();
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
    backup.clear();
    backup.add_bool("continue", false);
    save_block_to_file("backup.blk", backup);
  }
}