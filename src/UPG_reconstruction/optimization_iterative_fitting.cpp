#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>
#include "sdf_node.h"

namespace upg
{
  // creature is a unit of genetic algorithm
  // it is a set of parameters P with some additional info and
  // values of loss(P) and fitness function
  struct Creature
  {
    Creature() = default;
    Creature(const OptParams &p, int epoch)
    {
      params = p;
      local_opt_num = 0;
    }
    Creature(const OptParams &p, float l, int epoch)
    {
      params = p;
      loss = l;
      local_opt_num = 0;
    }
    OptParams params;
    std::shared_ptr<UPGOptimizer> local_optimizer;
    float loss = 1e9;
    float fitness = -1;
    int local_opt_num = 0;
  };

  enum InitStrategy
  {
    UNIFORM,
    LOG_UNIFORM
  };

  Creature initialize_creature(const std::vector<glm::vec2> &borders, const std::vector<int> &position_indices,
                               InitStrategy init_strategy = InitStrategy::LOG_UNIFORM)
  {
    Creature c;
    c.params.resize(borders.size());

    // initializing all parameters with some initialization strategy
    // position parameters are initialized too, but they'll be rewritten later
    for (int j = 0; j < borders.size(); j++)
    {
      if (init_strategy == InitStrategy::LOG_UNIFORM)
      {
        float t1 = log(borders[j].x);
        float t2 = log(borders[j].y);
        float l = t2 - t1;
        // logerr("%f %f t1 t2 %f %f", borders[j].x, borders[j].y, t1 , t2);
        c.params[j] = exp(urand(t1 + 0.25 * l, t2 - 0.25 * l));
      }
      else
        c.params[j] = CLAMP(borders[j].x + urand() * (borders[j].y - borders[j].x), borders[j].x, borders[j].y);
    }

    // special initialization for position parameters
    // We are trying to spawn primitives far from each other

    int primitives_cnt = position_indices.size();
    std::vector<glm::vec3> positions(primitives_cnt);
    positions[0] = glm::vec3(urand(borders[position_indices[0] + 0].x, borders[position_indices[0] + 0].y),
                             urand(borders[position_indices[0] + 1].x, borders[position_indices[0] + 1].y),
                             urand(borders[position_indices[0] + 2].x, borders[position_indices[0] + 2].y));
    for (int i = 1; i < primitives_cnt; i++)
    {
      int tries = 100;
      float best_dist = 0;
      glm::vec3 best_pos = positions[0];
      for (int j = 0; j < tries; j++)
      {
        float d = 1e9;
        glm::vec3 p = glm::vec3(urand(borders[position_indices[i] + 0].x, borders[position_indices[i] + 0].y),
                                urand(borders[position_indices[i] + 1].x, borders[position_indices[i] + 1].y),
                                urand(borders[position_indices[i] + 2].x, borders[position_indices[i] + 2].y));
        for (int k = 0; k < i; k++)
          d = std::min(d, length(p - positions[k]));
        if (d > best_dist)
        {
          best_dist = d;
          best_pos = p;
        }
      }
      positions[i] = best_pos;
    }

    for (int i = 0; i < primitives_cnt; i++)
    {
      c.params[position_indices[i] + 0] = positions[i].x;
      c.params[position_indices[i] + 1] = positions[i].y;
      c.params[position_indices[i] + 2] = positions[i].z;
    }
    return c;
  }

  //simple mutation
  void mutation(const std::vector<glm::vec2> &borders, OptParams &params, float mutation_chance, float mutation_power)
  {
    Normal normal_gen = Normal(0, mutation_power);

    for (int i = 0; i < MAX(1, mutation_chance * params.size()); i++)
    {
      /*
      int id = urandi(0, params.size());
      std::vector<float> mutation_power_chances(stat[0].size(), 0);
      float bin_f = stat[0].size() * ((params[id] - borders[id].x) / (borders[id].y - borders[id].x));

      for (int j = 0; j < stat[0].size(); j++)
      {
        mutation_power_chances[j] = (j == 0) ? 0 : mutation_power_chances[j - 1];
        mutation_power_chances[j] += (1 / (stat[id][j].second + 1e-6)) * normal_pdf((float)j / (stat[0].size()), bin_f, mutation_power);
      }

      
      float rnd = urand(0, mutation_power_chances.back());
      int mut_id = 0;
      while (rnd > mutation_power_chances[mut_id])
        mut_id++;
      mut_id += urand();
      float t = borders[id].x + ((float)mut_id/(stat[0].size()))*(borders[id].y - borders[id].x);
      t = CLAMP(t, borders[id].x, borders[id].y);
      */    
      int id = urandi(0, params.size());
      float t = borders[id].y + 1;
      while (t >= borders[id].y || t <= borders[id].x)
        t = CLAMP(params[id] + normal_gen.get() * (borders[id].y - borders[id].x), borders[id].x, borders[id].y);

      params[id] = t;
    }
  }

  //one-dot crossover
  OptParams crossover(const OptParams &p1, const OptParams &p2)
  {
    int id = urandi(0, p1.size());
    OptParams p = p1;
    for (int i = id; i < p1.size(); i++)
      p[i] = p2[i];
    return p;
  }

  class UPGIterativeFittingOptimizer : public UPGOptimizer
  {
  private:
    void initialize_population()
    {
      population.clear();
      best_population.clear();
      population.resize(population_size);
      best_population.resize(best_params_size);

      for (int i=0;i<population_size;i++)
      {
        population[i] = initialize_creature(borders, position_indices);
        evaluate(population[i]);
      }
    }

    float normal_pdf(float x, float x0, float sigma)
    {
      return (1/sqrt(2*PI*sigma*sigma)) * exp(-SQR(x-x0)/(2*sigma*sigma));
    }

    void local_search(Creature &c, int iterations, bool precise)
    {
      if (!c.local_optimizer || precise)
      {
        UPGReconstructionResult res;
        res.structure = structure;
        res.parameters.p = func->opt_params_to_gen_params(c.params, pd);

        c.local_optimizer = get_optimizer_adam(func, precise ? precise_opt_block : local_opt_block, res);
      }
      c.local_optimizer->optimize(iterations);
      UPGReconstructionResult lo_res = c.local_optimizer->get_best_results()[0];

      c.params = func->gen_params_to_opt_params(lo_res.parameters.p, pd);
      c.loss = lo_res.loss_optimizer;
      c.local_opt_num++;
      diff_function_calls += iterations;
      register_new_param(c);
    }

    void register_new_param(const Creature &c)
    {
      float worst_best_res = best_population[0].loss;
      int worst_index = 0;
      for (int i=1; i<best_params_size; i++)
      {
        if (best_population[i].loss > worst_best_res)
        {
          worst_index = i;
          worst_best_res = best_population[i].loss;
        }
      }
      
      if (worst_best_res > c.loss)
      {
        if (verbose)
        {
          debug("%d new best (%f) [", no_diff_function_calls, c.loss);
          for (int i=0;i<c.params.size();i++)
            debug("%f ", c.params[i]);
          debug("]\n");
        }
        best_population[worst_index] = c;
      }

      result_bins[CLAMP((int)(result_bin_count*c.loss), 0, result_bin_count-1)]++;
      for (int i=0;i<borders.size();i++)
      {
        int bin = stat[i].size()*((c.params[i] - borders[i].x)/(borders[i].y-borders[i].x));
        bin = CLAMP(bin, 0, stat[i].size()-1);
        stat[i][bin].first++;
        stat[i][bin].second = std::min(stat[i][bin].second, (double)(c.loss));
      }
    }

    void evaluate(Creature &c)
    {
      c.loss = func->f_no_grad(gen.get(), pd, c.params);
      register_new_param(c);
      no_diff_function_calls++;
    }
  
    std::pair<int, int> choose_parents_tournament(int tournament_size)
    {
      std::vector<int> indexes = std::vector<int>(tournament_size, 0);
      for (int i=0;i<tournament_size;i++)
        indexes[i] = urandi(0, population_size);
      std::sort(indexes.begin(), indexes.end(), [this](const int & a, const int & b) -> bool{    
            return population[a].loss < population[b].loss;});

      return std::pair<int, int>(indexes[0], indexes[1]);
    }

    float pos_fitness_function(int pos)
    {
      return 1.0/(pos+1);
    }

    void sort_and_calculate_fitness()
    {
      std::sort(population.begin(), population.end(), [](const Creature & a, const Creature & b) -> bool{    
                return a.loss < b.loss;});
      for (int i=0;i<population_size;i++)
        population[i].fitness = pos_fitness_function(i);
    }

    void print_result_bins(int *bins)
    {
      int total_cnt = 0;
      for (int i=0;i<result_bin_count;i++)
        total_cnt+=bins[i];
      for (int i=0;i<10;i++)
      {
        int cnt = 0;
        for (int j=0;j<1;j++)
          cnt += bins[i+j];
        debug("[%.4f-%.4f] %d (%.2f %)\n",(float)i/result_bin_count, (float)(i+1)/result_bin_count, cnt, 100.f*cnt/total_cnt);
      }
      for (int i=10;i<100;i+=10)
      {
        int cnt = 0;
        for (int j=0;j<10;j++)
          cnt += bins[i+j];
        debug("[%.4f-%.4f] %d (%.2f %)\n",(float)i/result_bin_count, (float)(i+10)/result_bin_count, cnt, 100.f*cnt/total_cnt);
      }
      for (int i=100;i<result_bin_count;i+=100)
      {
        int cnt = 0;
        for (int j=0;j<100;j++)
          cnt += bins[i+j];
        debug("[%.4f-%.4f] %d (%.2f %)\n",(float)i/result_bin_count, (float)(i+100)/result_bin_count, cnt, 100.f*cnt/total_cnt);
      }
    }

  public:
    UPGIterativeFittingOptimizer(UPGOptimizableFunction *_func, 
                                 const Block &settings, const UPGStructure &_structure) :
    UPGOptimizer(_func)
    {
      structure = _structure;
      gen = func->get_generator(structure);
      pd = func->get_full_parameters_description(gen.get());

      verbose = settings.get_bool("verbose");
      finish_thr = settings.get_double("finish_threshold");
      budget = settings.get_int("iterations", budget);
      local_opt_block.set_bool("verbose", verbose);
      local_opt_block.set_bool("save_intermediate_images", false);
      local_opt_block.set_double("learning_rate", local_learning_rate);
      precise_opt_block.set_bool("verbose", verbose);
      precise_opt_block.set_bool("save_intermediate_images", false);
      precise_opt_block.set_double("learning_rate", 0.02);
    }

    UPGIterativeFittingOptimizer(UPGOptimizableFunction *_func, const Block &settings,
                                 const UPGStructure &_structure, const UPGParametersRaw &_params) :
    UPGIterativeFittingOptimizer(_func, settings, _structure)
    {
      ref_params = _params.p;
      is_test_stand = true;
    }

    virtual void optimize(int iters = -1) override
    {
      assert(dynamic_cast<SdfGenInstance*>(gen.get())); //currently not working with meshes

      if (is_test_stand)
        optimization_test_stand(iters);
      else
        optimize_part_based(iters);
    }

    void set_parameters_info(const ParametersDescription &param_desc, const std::vector<UPGPart> &parts)
    {
      borders.clear();
      position_indices.clear();
      opt_to_gen_indices.clear();
      gen_to_opt_indices.clear();

      int gi = 0;
      int oi = 0;
      for (const auto &p : pd.get_block_params())
      {
        for (const auto &param_info : p.second.p)
        {
          if (param_info.type != ParameterType::CONST)
          {
            borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
            gen_to_opt_indices.push_back(oi);
            opt_to_gen_indices.push_back(gi);
            oi++;
          }
          else
            gen_to_opt_indices.push_back(-1);
          gi++;
        }
      }

      for (auto &p : parts)
      {
        if (p.position_index >= 0 && gen_to_opt_indices[p.position_index] >= 0)
        {
          position_indices.push_back(gen_to_opt_indices[p.position_index]);
          logerr("added PI %d", position_indices.back());
        }
      }
    }

    void optimize_part_based(int iters)
    {
      population_size = 1000;
      int step_budget = 25000;
      float quality_thr = 0.001;
      
      iters = iters > 0 ? iters : budget;
      int steps = ceil(iters/step_budget);
      std::vector<UPGPart> parts = get_sdf_parts(structure);
      std::vector<bool> fixed_parts(parts.size(), false);
      std::vector<bool> fixed_params(pd.get_total_params_count(), false);
      int parts_left = parts.size();
      int step = 0;

      //pd.print_info();

      while (step < steps && parts_left > 0)
      {
        //set actual info about parameters and parts
        set_parameters_info(pd, parts);

        //optimize current setup using memetic optimizer
        optimize_memetic(step_budget);
        local_search(best_population[0], 300, true);

        //check all part to find ones that should be fixed
        auto best_res = get_best_results()[0];
        debug("fixed [ ");
        parts_left = parts.size();
        for (int i=0;i<parts.size();i++)
        {
          float q = func->estimate_positioning_quality(structure, parts[i], best_res.parameters.p, 0.01, 100);
          if (q > quality_thr)
          {
            fixed_parts[i] = true;
            parts_left--;
          }
          debug("%d ",(int)fixed_parts[i]);
        }
        debug("] best %.8f\n", best_population[0].loss);

        //make parameters in fixed parts constant
        ParametersDescription old_pd = pd;
        for (int i=0;i<parts.size();i++)
        {
          if (!fixed_parts[i])
            continue;
          int p_id = parts[i].p_range.first;
          for (unsigned n_id = parts[i].s_range.first; n_id < parts[i].s_range.second; n_id++)
          {
            assert(pd.get_block_params().find(n_id) != pd.get_block_params().end());
            for (auto &param : pd.get_block_params().at(n_id).p)
            {
              //logerr("fixed param %d %d %s", n_id, p_id, param.name.c_str());
              assert(p_id < parts[i].p_range.second);
              param.type = ParameterType::CONST;
              param.value = best_res.parameters.p[p_id];
              fixed_params[p_id] = true;
              p_id++;
            }
          }
        }

        //set constant values for these fixed parameters for the whole population
        for (auto &c : population)
        {
          //population was coded with old parameters description (less constants, more opt parameters)
          auto gen_params = func->opt_params_to_gen_params(c.params, old_pd); 
          c = Creature(func->gen_params_to_opt_params(gen_params, pd), 0);
        }

        step++;
        if (best_population[0].loss < finish_thr)
          break;
      }
    }

    void optimize_memetic(int iters)
    {
      borders.clear();

      for (const auto &p : pd.get_block_params())
      {
        for (const auto &param_info : p.second.p)
        {
          if (param_info.type != ParameterType::CONST)
          {
            borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
          }
        }
      }
  
      stat.resize(borders.size());
      for (auto &a : stat)
        for (auto &p : a)
          p = {0, 1.0};
      if (iters > 0)
        budget = iters;
      epoch = 0;
      no_diff_function_calls = 0;
      diff_function_calls = 0;

      if (population.empty())
        initialize_population();

      while (no_diff_function_calls + diff_function_calls < budget)
      {
        int elites_count = elites_fraction*population_size;
        std::vector<Creature> new_population(population_size);

        sort_and_calculate_fitness();
        for (int i=0;i<elites_count;i++)
          new_population[i] = population[i];

        for (int i=elites_count;i<population_size;i++)
        {
          float time_q = (no_diff_function_calls + diff_function_calls + 0.0) / budget;
          float mutation_chance = (1-time_q)*max_mutation_chance + time_q*min_mutation_chance;
          float mutation_power = (1-time_q)*max_mutation_power + time_q*min_mutation_power;
          std::vector<float> chances = {0.1, 0.1 + 0.5, 0.1 + 0.5 + 0.4};
          float action_rnd = urand(0, chances.back());
          if (action_rnd < chances[0])
          {
            //create new random creature from scratch
            new_population[i] = initialize_creature(borders, position_indices);
            evaluate(new_population[i]);
          }
          else if (action_rnd < chances[1])
          {
            OptParams p = population[urandi(0, elites_count)].params;
            mutation(borders, p, mutation_chance, urand()*mutation_power);
            new_population[i] = Creature(p, epoch);
            evaluate(new_population[i]);
          }
          else
          {
            auto parents = choose_parents_tournament(tournament_size);
            OptParams p = crossover(population[parents.first].params, population[parents.second].params);
            mutation(borders, p, mutation_chance, urand()*mutation_power);
            new_population[i] = Creature(p, epoch);
            evaluate(new_population[i]);
          }
        }

        for (int i=0; i<result_bin_count; i++)
          current_bins[i] = 0;
        for (auto &c : new_population)
          current_bins[CLAMP((int)(result_bin_count*c.loss), 0, result_bin_count-1)]++;
        population = new_population;

        if (verbose)
        {
          debug("EPOCH %d (%d+%d/%d)\n", epoch, no_diff_function_calls, diff_function_calls, budget);
          for (auto &c : best_population)
          {
            debug("best res %.6f [", c.loss);
            for (int i=0;i<c.params.size();i++)
              debug("%f ", c.params[i]);
            debug("]\n");
          }
          debug("\n");
          /*
          for (auto &arr : stat)
          {
            debug("param stat: ");
            for (auto &p : arr)
              debug("[%d %f]",p.first, (float)p.second);
            debug("\n");
          }
          */
        }

        sort_and_calculate_fitness();
        int good_soulutions = 0;
        for (int i=0;i<population_size;i++)
        {
          if (population[i].loss < good_soulution_thr)
            good_soulutions++;
          else
            break;
        }
        std::vector<float> chances(good_soulutions,0);
        chances[0] = 1;
        for (int i=1;i<good_soulutions;i++)
          chances[i] = chances[i-1] + 1.0/(i+1);
        std::vector<int> indices_to_opt(local_opt_count, 0);
        for (int i=0;i<local_opt_count;i++)
        {
          float rnd = urand(0, chances.back());
          for (int j=0;j<good_soulutions;j++)
          {
            if (rnd<chances[j])
            {
              indices_to_opt[i] = j;
              break;
            }
          }
        }

        for (auto &i : indices_to_opt)
          local_search(population[i], local_opt_iters, false);

        epoch++;
        if (best_population[0].loss < finish_thr)
          break;
      }
    }
    virtual std::vector<UPGReconstructionResult> get_best_results() override
    {
      UPGReconstructionResult res;
      res.structure = structure;
      {
        res.parameters.p = func->opt_params_to_gen_params(best_population[0].params, pd);
        res.loss_optimizer = best_population[0].loss;
      }
      return {res};
    }
  private:
    //settings
    int population_size = 5000;
    int best_params_size = 1;
    int budget = 200'000; //total number of function calls
    float max_mutation_chance = 0.5;
    float min_mutation_chance = 0.1;
    float max_mutation_power = 0.5;
    float min_mutation_power = 0.3;
    int tournament_size = 200;
    int local_opt_count = 50;
    int local_opt_iters = 50;
    float local_learning_rate = 0.05;
    float good_soulution_thr = 0.001;
    float elites_fraction = 0.01;
    float finish_thr = 0;
    bool verbose = false;
    bool is_test_stand = false;

    //info about structure of parameters
    std::vector<glm::vec2> borders;      //size equals total size of OptParams vector
    std::vector<int> position_indices;   //indices in OptVector that represents positions of primitives (one for every non-fixed primitive)
    std::vector<int> gen_to_opt_indices; //-1 mean that this parameter is fixed and shouldn't be optimized
    std::vector<int> opt_to_gen_indices; //all elements >= 0

    //generator-specific stuff
    UPGStructure structure;
    std::shared_ptr<UniversalGenInstance> gen;
    ParametersDescription pd;
    Block local_opt_block, precise_opt_block;
    std::vector<std::array<std::pair<int, double>, 64>> stat;
    std::vector<float> ref_params;

    //GA state
    std::vector<Creature> population;
    std::vector<Creature> best_population;

    int epoch = 0;
    int no_diff_function_calls = 0;
    int diff_function_calls = 0;

    //statistics
    constexpr static int result_bin_count = 1000;
    int result_bins[result_bin_count] = {};
    int current_bins[result_bin_count] = {};
  
    class Oracle
    {
      std::vector<float> reference_params;
      std::vector<glm::vec2> borders;
    public:
      Oracle() = default;
      Oracle(std::vector<float> &_reference_params, const std::vector<glm::vec2> &_borders)
      {
        assert(_reference_params.size() == _borders.size());
        reference_params = _reference_params;
        borders = _borders;
      }
      Creature initialize_creature(float hit_chance, float variance)
      {
        Normal normal(0, variance);
        Creature c;
        c.params.resize(borders.size());
        for (int j=0;j<borders.size();j++)
        {
          if (urand() < hit_chance)
          {
            c.params[j] = CLAMP(reference_params[j] + normal.get()*(borders[j].y - borders[j].x), borders[j].x, borders[j].y);
          }
          else
          {
            if (j/3%2==1)
            {
              float t1 = log(borders[j].x);
              float t2 = log(borders[j].y);
              float l = t2-t1;
              //logerr("%f %f t1 t2 %f %f", borders[j].x, borders[j].y, t1 , t2);
              c.params[j] = exp(urand(t1+0.25*l, t2-0.25*l));
            }
            else
              c.params[j] = CLAMP(borders[j].x + urand()*(borders[j].y - borders[j].x), borders[j].x, borders[j].y);
          }
        }

        int primitives_cnt = borders.size()/6;
        std::vector<glm::vec3> positions(primitives_cnt);
        positions[0] = glm::vec3(urand(borders[0].x, borders[0].y), urand(borders[1].x, borders[1].y), urand(borders[2].x, borders[2].y));
        for (int i=1;i<primitives_cnt;i++)
        {
          int tries = 100;
          float best_dist = 0;
          glm::vec3 best_pos = positions[0];
          for (int j=0;j<tries;j++)
          {
            float d = 1e9;
            glm::vec3 p = glm::vec3(urand(borders[0].x, borders[0].y), urand(borders[1].x, borders[1].y), urand(borders[2].x, borders[2].y));
            for (int k=0;k<i;k++)
              d = std::min(d, length(p-positions[k]));
            if (d > best_dist)
            {
              best_dist = d;
              best_pos = p;
            }
          }
          positions[i] = best_pos;
        }

        for (int j=0;j<primitives_cnt;j++)
        {
          c.params[6*j] = positions[j].x;
          c.params[6*j+1] = positions[j].y;
          c.params[6*j+2] = positions[j].z;
        }

        return c;
      }
    };

    void optimization_test_stand(int iters)
    {
      borders.clear();

      for (const auto &p : pd.get_block_params())
      {
        for (const auto &param_info : p.second.p)
        {
          if (param_info.type != ParameterType::CONST)
          {
            borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
          }
        }
      }

      stat.resize(borders.size());
      for (auto &a : stat)
        for (auto &p : a)
          p = {0, 1.0};
      if (iters > 0)
        budget = iters;
      epoch = 0;
      no_diff_function_calls = 0;
      diff_function_calls = 0;
      best_population.resize(best_params_size);

      Oracle oracle(ref_params, borders);

      int k = 0;
      int found = 0;

      // #pragma omp parallel for
      for (int i = 0; i < 1000; i++)
      {
        Creature c = oracle.initialize_creature(0, 0.5);
        local_search(c, 200, true);
        if (c.loss <= 1e-6)
          found++;
        // logerr("%d res = %f. best = %f", k, c.loss, best_population[0].loss);
        k++;
      }
      logerr("found %d/%d", found, k);
      return;
    }

  };

  std::shared_ptr<UPGOptimizer> get_optimizer_iterative_fitting(UPGOptimizableFunction *_func, 
                                                                const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGIterativeFittingOptimizer>(_func, settings, _structure);
  }
  std::shared_ptr<UPGOptimizer> optimization_test_stand(UPGOptimizableFunction *_func, const Block &settings, 
                                                        const UPGStructure &_structure, const UPGParametersRaw &_params)
  {
    return std::make_shared<UPGIterativeFittingOptimizer>(_func, settings, _structure, _params);
  }
}