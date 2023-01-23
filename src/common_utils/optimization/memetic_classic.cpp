#include "optimization.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"
#include <algorithm>

namespace opt
{
  struct Solution
  {
    std::vector<float> params = {};
    float value = 1000;
    float qa_value = 1000;
    int local_improvements = 0;
  };

  void MemeticClassic::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                                opt_func_vector &F_reg)
  {
    float mutation_chance = settings.get_double("mutation_chance", 0.2);
    float mutation_power = settings.get_double("mutation_power", 0.3);
    int tournament_size = settings.get_int("tournament_size", 16);
    int population_size = settings.get_int("population_size", 32);
    int total_function_calls = settings.get_int("total_function_calls", 4800);
    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;
    int local_search_iterations = settings.get_int("local_search_iterations", 50);
    float local_search_learning_rate = settings.get_double("local_search_learning_rate", 0.01);
    int start_search_iterations = settings.get_int("start_search_iterations", 25);
    float start_search_learning_rate = settings.get_double("start_search_learning_rate", 0.05);
    float recreation_diversity_thr = settings.get_double("recreation_diversity_thr", 1.75);
    float depth_reg_q = settings.get_double("depth_reg_q", 0);
    int budget = 0;

    auto mutate = [&](const std::vector<float> &base) -> std::vector<float>
    {
      std::vector<float> res = base;
      int id = urandi(0, base.size());
      res[id] = urand(min_X[id], max_X[id]);
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

    auto choose_parents_tournament = [&](std::vector<Solution> &solutions) -> std::pair<int, int>
    {
      std::vector<int> indexes = std::vector<int>(tournament_size, 0);
      for (int i=0;i<tournament_size;i++)
        indexes[i] = urandi(0, solutions.size());
      std::sort(indexes.begin(), indexes.end(), [&](const int & a, const int & b) -> bool{    
            return solutions[a].value < solutions[b].value;});

      return std::pair<int, int>(indexes[0], indexes[1]);
    };

    auto initialize_population = [&](std::vector<Solution> &solutions)
    {
      solutions = std::vector<Solution>(population_size, {std::vector<float>(min_X.size(), 0), 1000, 1000, 0});
      for (int i=0;i<population_size;i++)
      {
        for (int j=0;j<min_X.size();j++)
          solutions[i].params[j] = urand(min_X[j], max_X[j]);
      }

      std::vector<std::vector<float>> presets_1 = {
        {3    , 3    , 3    , 3    , 3    , 3    , 3    , 3    , 3    , 0.75 , 1},
        {3    , 3    , 3    , 3    , 3    , 3    , 3    , 3    , 3    , 1    , 1},
        {4    , 4    , 4    , 4    , 4    , 4    , 4    , 4    , 4    , 0.75 , 1},
        {4    , 4    , 4    , 4    , 4    , 4    , 4    , 4    , 4    , 1    , 1},
        {1.879, 1.888, 2.867, 3.070, 3.333, 3.533, 3.746, 3.899, 3.92 , 0.690, 0},
        {5.5  , 5.8  , 6.1  , 6.4  , 6.7  , 7    , 7.3  , 7.6  , 7.9  , 0.1  , 0}
      };

      std::vector<std::vector<float>> presets_2 = {
        {0.05, 0.5, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33},
        {0.05, 0.6, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33},
        {0.05, 0.4, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33},
        {0.05, 0.5, 0.2 , 0.2 , 0.2 , 0.2 , 0.2 , 0.2 , 0.2 },
      };

      std::vector<std::vector<float>> presets_3 = {
        {0.1, 0, 0.004, 0.3, 3.141, -0.250},
        {0.1, 0, 0.004, 0.0, 0.0, 0.1},
        {0.1, 0.2, -0.2, 0.05, 1, 0.0},
        {0.1, 0.2, -0.2, 0.1, 2, 0.3},
        {0.1, 0.2, -0.2, 0.2, 0.0, 0.03},
        {-0.2, 0.60, 0.2, 0.4, 3.14, 0.0},
        {-0.2, 0.5, 0.5, 0.4, 3.14, 0.0},
        {-0.2, 0.60, 0.2, 0.3, 0, 0.0},
        {-0.2, 0.60, 0.2, 0.4, 3.14, 0.0},
        {0, 0.5, 0.5, 0.4, 0, 0.0},
        {0, 0.60, 0.2, 0.3, 0, 0.0}
      };

      std::vector<std::vector<float>> all_presets;
      for (auto &p1 : presets_1)
      {
        for (auto &p2 : presets_2)
        {
          for (auto &p3 : presets_3)
          {
            all_presets.emplace_back();
            for (auto &v1 : p1)
              all_presets.back().push_back(v1);
            for (auto &v2 : p2)
              all_presets.back().push_back(v2);
            for (auto &v3 : p3)
              all_presets.back().push_back(v3);
          }
        }
      }

      for (int i=0; i<population_size; i++)
      {
        solutions[i].params = all_presets[(int)urandi(0, all_presets.size())];
      }
    };

    auto local_search = [&](std::vector<float> start_params, int iterations, float lr) -> std::pair<float, std::vector<float>>
    {
      Optimizer *opt = new Adam();
      Block adam_settings;
      adam_settings.add_arr("initial_params", start_params);
      adam_settings.add_double("learning_rate", lr);
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
        debug("%d new best[1] %.4f\n", budget, best_result);
      }
      delete opt;
      return res;
    };

    auto calc_diversity = [&](std::vector<Solution> &population) -> float
    {
      float diversity = 0;
      int cnt = 0;
      for (int i1=0;i1<population.size();i1++)
      {
        for (int i2=i1+1;i2<population.size();i2++)
        {
          float d = 0;
          for (int j=0;j<min_X.size();j++)
            d += SQR((population[i1].params[j] - population[i2].params[j])/(max_X[j] - min_X[j]));
          diversity += sqrtf(d) / min_X.size();
          cnt++;
        }
      }
      diversity /= cnt;
      return diversity;
    };

    auto load_solutions = [](Block &backup, std::vector<Solution> &solutions)
    {
      solutions.clear();
      for (int i=0;i<backup.size();i++)
      {
        Block *s_block = backup.get_block(i);
        if (s_block)
        {
          solutions.emplace_back();
          s_block->get_arr("params", solutions.back().params);
          solutions.back().value = s_block->get_double("value");
          solutions.back().qa_value = s_block->get_double("qa_value");
          solutions.back().local_improvements = s_block->get_int("local_improvements");
        }
      }
    };

    auto save_solutions = [](Block &backup, std::vector<Solution> &solutions)
    {
      for (int i=0;i<solutions.size();i++)
      {
        Block s_block;
        s_block.set_arr("params", solutions[i].params);
        s_block.set_double("value", solutions[i].value);
        s_block.set_double("qa_value", solutions[i].qa_value);
        s_block.set_int("local_improvements", solutions[i].local_improvements);

        backup.add_block("solution_"+std::to_string(i), &s_block);
      }
    };
    std::vector<Solution> solutions;

    Block backup;
    if (load_block_from_file("backup.blk", backup) && backup.get_bool("continue"))
    {
      budget = backup.get_int("budget");
      load_solutions(backup, solutions);
      for (auto &s : solutions)
      {
        if (s.value < best_result)
        {
          best_result = s.value;
          best_params = s.params;
          debug("%d new best[0] %.4f\n", budget, best_result);
        }
      }
    }
    else
    {
      initialize_population(solutions);
      for (int i = 0; i < solutions.size(); i++)
      {
        auto res = local_search(solutions[i].params, start_search_iterations, start_search_learning_rate);

        solutions[i].params = res.second;
        solutions[i].value = res.first;
        solutions[i].qa_value = 1e-3 / (res.first * res.first);
        solutions[i].local_improvements = 1;
      }
    }
    float initial_diversity = calc_diversity(solutions);
    debug("initial diversity %.3f\n", initial_diversity);

    while (budget < total_function_calls)
    {
      backup.clear();
      backup.add_bool("continue", true);
      backup.add_int("budget", budget);
      save_solutions(backup, solutions);
      save_block_to_file("backup.blk", backup);

      std::pair<float, std::vector<float>> res;
      if (urand() < 0.2*((float)budget)/total_function_calls)
      {
        //crossover + mutation
        auto parents = choose_parents_tournament(solutions);
        std::vector<float> new_solution = one_dot_crossover(solutions[parents.first].params, solutions[parents.second].params);
        new_solution = mutate(new_solution);
        res = local_search(new_solution, local_search_iterations, local_search_learning_rate);
      }
      else
      {
        int id = urandi(0, solutions.size());
        std::vector<std::vector<float>> new_solution = {solutions[id].params};
        new_solution[0] = mutate(new_solution[0]);
        auto fres = F(new_solution);
        res.first = fres[0].first;
        res.second = new_solution[0];
        if (res.first < best_result)
        {
          best_result = res.first;
          best_params = res.second;
          debug("%d new best[2] %.4f\n", budget, best_result);
        }
        budget++;
      }
      int worst_idx = -1;
      float worst_val = res.first;
      for (int i=0;i<solutions.size();i++)
      {
        if (solutions[i].value > worst_val)
        {
          worst_val = solutions[i].value;
          worst_idx = i;

        }
      }
      if (worst_idx >= 0)
      {
        solutions[worst_idx] = Solution{res.second, res.first, 1e-3f / (res.first * res.first), 1};
      }

      //further local improvement
      if (urand() < ((float)budget)/total_function_calls)
      {
        float qav_sum = 0;
        for (auto &s : solutions)
          qav_sum += s.qa_value;

        float rnd = urand(qav_sum);
        int improve_idx = 0;
        for (int i=0;i<solutions.size();i++)
        {
          if (rnd < solutions[i].qa_value)
          {
            improve_idx = i;
            break;
          }
          else
          {
            rnd -= solutions[i].qa_value;
          }
        }
        auto res = local_search(solutions[improve_idx].params, local_search_iterations, local_search_learning_rate);
        if (solutions[improve_idx].value > res.first)
        {
          logerr("improving existing solution (%d times) %.4f --> %.4f", solutions[improve_idx].local_improvements, solutions[improve_idx].value,
                 res.first);
          solutions[improve_idx].params = res.second;
          solutions[improve_idx].value = res.first;
          solutions[improve_idx].qa_value = 1e-3 / (res.first * res.first);
          solutions[improve_idx].local_improvements = 1;
        }
        else
        {
          //if we improve solution, there is unlikely that it will be further improved
          solutions[improve_idx].local_improvements++;
          float base_qa = 1e-3 /(solutions[improve_idx].value * solutions[improve_idx].value);
          solutions[improve_idx].qa_value = base_qa / solutions[improve_idx].local_improvements;
        }
      }

      //diversity calculation and population restart
      float diversity = calc_diversity(solutions);
      if (initial_diversity/(diversity + 1e-6) > recreation_diversity_thr)
      {
        //recreate population and place best solution in it
        debug("recreating population, diversity = %.2f\n", diversity);
        initialize_population(solutions);
        solutions[0] = Solution{best_params, best_result, 1e-3f / (best_result * best_result), 1};
        for (int i=0; i<solutions.size(); i++)
        {
          auto res = local_search(solutions[i].params, start_search_iterations, start_search_learning_rate);
          solutions[i] = Solution{res.second, res.first, 1e-3f / (res.first * res.first), 1};
        }
        initial_diversity = calc_diversity(solutions);  
        debug("initial diversity %.3f\n", initial_diversity);      
      }
      
    }
    backup.clear();
    backup.add_bool("continue", false);
    save_block_to_file("backup.blk", backup);
    remove("config/backup.blk");
  }
}