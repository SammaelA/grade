#include "optimization.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <bits/stdc++.h>

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
                                init_params_func &get_init_params)
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
    int min_iter_between_recreations = settings.get_int("min_iter_between_recreations", 250);
    
    //sometimes we want to fix specific parameters during memetic optimization
    //Main reason for it is for multi-stage optimization when we want to preserve
    //values, already optimized on previous stage
    //To fix the parameters we need their values, that's why initial_params appear here

    std::vector<float> derivatives_mult;
    settings.get_arr("derivatives_mult", derivatives_mult);
    if (derivatives_mult.size()!=min_X.size())
      derivatives_mult = std::vector<float>(min_X.size(), 1);
    std::vector<float> initial_params;
    settings.get_arr("initial_params", initial_params);
    if (initial_params.empty())
    {
      initial_params = std::vector<float>(min_X.size(), 0);
      for (int i = 0; i < initial_params.size(); i++)
      {
        initial_params[i] = min_X[i] + ((urand(0, 1) + urand(0, 1)) / 2) * (max_X[i] - min_X[i]);
      }
    }

    int budget = 0;
    Normal normal_gen = Normal(0, mutation_power);

    auto initialize = [&]() -> std::vector<float>
    {
      std::vector<float> p = get_init_params();
      for (int n=0;n<min_X.size();n++)
      {
        if (derivatives_mult[n] <= 0)
          p[n] = initial_params[n];
      }
      return p;
    };

    auto mutate = [&](const std::vector<float> &base) -> std::vector<float>
    {
      std::vector<float> res = base;
      for (int i=0;i<MAX(1, mutation_chance * res.size());i++)
      {
        int id = urandi(0, base.size());
        if (derivatives_mult[id] > 0)
        {
          float t = max_X[id] + 1;
          while (t >= max_X[id] || t <= min_X[id])
            t = CLAMP(res[id] + normal_gen.get()*(max_X[id] - min_X[id]), min_X[id], max_X[id]);
          res[id] = t;
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

    auto choose_parents_tournament = [&](std::vector<Solution> &solutions) -> std::pair<int, int>
    {
      std::vector<int> indexes = std::vector<int>(tournament_size, 0);
      for (int i=0;i<tournament_size;i++)
        indexes[i] = urandi(0, solutions.size());
      std::sort(indexes.begin(), indexes.end(), [&](const int & a, const int & b) -> bool{    
            return solutions[a].value < solutions[b].value;});

      return std::pair<int, int>(indexes[0], indexes[1]);
    };

    auto local_search = [&](std::vector<float> start_params, int iterations, float lr) -> std::pair<float, std::vector<float>>
    {
      Optimizer *opt = new Adam();
      Block adam_settings;
      adam_settings.add_arr("initial_params", start_params);
      adam_settings.add_double("learning_rate", lr);
      adam_settings.add_int("iterations", iterations);
      adam_settings.add_bool("verbose", false);
      adam_settings.add_arr("derivatives_mult", derivatives_mult);
      opt->optimize(F, min_X, max_X, adam_settings);
      budget += iterations;

      std::pair<float, std::vector<float>> res;
      res.second = opt->get_best_result(&(res.first));

      delete opt;
      return res;
    };

    auto initialize_population = [&](std::vector<Solution> &solutions, int st)
    {
      solutions = std::vector<Solution>(population_size, {std::vector<float>(min_X.size(), 0), 1000, 1000, 0});
      for (int i=0;i<population_size;i++)
        solutions[i].params = initialize();
      
        for (int i=st; i<solutions.size(); i++)
        {
          auto res = local_search(solutions[i].params, start_search_iterations, start_search_learning_rate);
          solutions[i] = Solution{res.second, res.first, 1e-3f / (res.first * res.first), 1};
          if (res.first < best_result)
          {
            best_result = res.first;
            best_params = res.second;
            debug("%d new best %.4f\n", budget, res.first);
          }
        }
      return;
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

        backup.set_block("solution_"+std::to_string(i), &s_block);
      }
    };

    auto block_hash = [](Block &b) -> uint64_t
    {
      std::hash<std::string> hash_string;
      std::string bs;
      save_block_to_string(bs, b);
      return hash_string(bs);
    };
    std::vector<Solution> solutions;
    Block backup;

    debug("Memetic classic: starting optimization\n");
    bool has_backup = load_block_from_file("backup.blk", backup) && backup.get_bool("continue");
    bool relevant_backup = false;
    if (has_backup)
    {
      relevant_backup = min_X.size() == backup.get_int("parameters_count") &&
                        block_hash(settings) == backup.get_uint64("settings_block_hash");
      if (!relevant_backup)
      {
        backup.clear();
        save_block_to_file("backup.blk", backup);
      }
    }
    if (relevant_backup)
    {
      debug("Continue from backup file\n");
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
      debug("Initializing population\n");
      initialize_population(solutions, 0);

      backup.set_bool("continue", true);
      backup.set_int("parameters_count", min_X.size());
      backup.set_uint64("settings_block_hash", block_hash(settings));
    }
    float initial_diversity = calc_diversity(solutions);
    debug("initial diversity %.3f\n", initial_diversity);

    std::vector<float> chances = {0.9, 0.7, 0.03, 0.00};
    std::vector<int> improvements = {0, 0, 0, 0};
    std::vector<int> best_stat = {0, 0, 0, 0};

    std::vector<float> best_value_history;
    int iterations_from_population_init = 0;

    std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
    int t_last_print = budget;
    int t_first_print = budget;

    while (budget < total_function_calls)
    {
      backup.set_int("budget", budget);
      save_solutions(backup, solutions);
      save_block_to_file("backup.blk", backup);

      int improvement_idx = 0;
      float rnd = urand();
      std::pair<float, std::vector<float>> res;

      if (rnd < chances[0])
      {
        //random mutation
        improvement_idx = 0;
        int id = urandi(0, solutions.size());
        std::vector<std::vector<float>> new_solution = {solutions[id].params};
        new_solution[0] = mutate(new_solution[0]);
        auto fres = F(new_solution);
        res.first = fres[0].first;
        res.second = new_solution[0];
        budget++;
      }
      else if (rnd < chances[0] + chances[1])
      {
        //crossover + mutation + local search
        improvement_idx = 1;
        auto parents = choose_parents_tournament(solutions);
        std::vector<float> new_solution = one_dot_crossover(solutions[parents.first].params, solutions[parents.second].params);
        new_solution = mutate(new_solution);
        res = local_search(new_solution, local_search_iterations, local_search_learning_rate);
      }
      else if (rnd < chances[0] + chances[1] + chances[2])
      {
        //local improvement
        improvement_idx = 2;
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
        res = local_search(solutions[improve_idx].params, local_search_iterations, local_search_learning_rate);
      }
      else
      {
        //new parameters set
        improvement_idx = 3;
        std::vector<float> new_params = initialize();
        res = local_search(new_params, start_search_iterations, start_search_learning_rate);
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
        improvements[improvement_idx]++;
      }

      if (res.first < best_result)
      {
        best_result = res.first;
        best_params = res.second;
        best_stat[improvement_idx]++;
        debug("%d new best - %.4f [%d %d %d %d]\n", budget, best_result, best_stat[0], best_stat[1], best_stat[2], best_stat[3]);
      }

      float diversity = calc_diversity(solutions);
      best_value_history.push_back(best_result);
      iterations_from_population_init++;

      bool diversity_dropped = initial_diversity/(diversity + 1e-6) > recreation_diversity_thr;
      bool quality_stall = best_value_history.size() > min_iter_between_recreations && 
                           best_value_history.back() - best_value_history[best_value_history.size() - min_iter_between_recreations] < 0.01;
      bool enough_budget_to_recreate = budget + population_size*start_search_iterations < total_function_calls;
      bool no_recent_recreation = iterations_from_population_init >= min_iter_between_recreations;
      //diversity calculation and population restart
      if ((diversity_dropped || quality_stall) &&
          enough_budget_to_recreate &&
          no_recent_recreation)
      {
        //recreate population and place best solution in it
        debug("recreating population, diversity = %.2f\n", diversity);
        iterations_from_population_init = 0;
        initialize_population(solutions, 1);
        solutions[0] = Solution{best_params, best_result, 1e-3f / (best_result * best_result), 1};
        initial_diversity = calc_diversity(solutions);  
        debug("initial diversity %.3f\n", initial_diversity);      
      }
      
      if (budget - t_last_print > 100)
      {
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        float ms_per_iter = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t_start).count()/(budget - t_first_print);
        int s = 1e-3*ms_per_iter*(total_function_calls - budget);
        debug("%d/%d best value %.4f diversity %.4f ETA: %d:%2d\n", budget, total_function_calls, best_result, diversity, s / 60, s % 60);
        t_last_print = budget;
      }
    }
    backup.clear();
    backup.add_bool("continue", false);
    save_block_to_file("backup.blk", backup);
    remove("config/backup.blk");
  }
}