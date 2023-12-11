#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>

namespace upg
{
  class UPGFixedStructureMemeticOptimizer : public UPGOptimizer
  {
  private:
    //creature is a unit of genetic algorithm
    //it is a set of parameters P with some additional info and 
    //values of loss(P) and fitness function
    struct Creature
    {
      Creature() = default;
      Creature(const OptParams &p, int epoch)
      {
        params = p;
        birth_epoch = epoch;
        local_opt_num = 0;
      }
      Creature(const OptParams &p, float l, int epoch)
      {
        params = p;
        loss = l;
        birth_epoch = epoch;
        local_opt_num = 0;
      }
      OptParams params;
      std::shared_ptr<UPGOptimizer> local_optimizer;
      float loss = 1e9;
      float fitness = -1;
      int birth_epoch = -1;
      int local_opt_num = 0;
    };
    //maps [0, +inf) to [0,1)
    float normalize_loss(float x)
    {
      return 2*(1/(1+expf(-CLAMP(x,0,10)))) - 1;
    }

    Creature initialize_creature()
    {
      Creature c;
      c.birth_epoch = epoch;
      c.params.resize(borders.size());
      for (int j=0;j<borders.size();j++)
        c.params[j] = borders[j].x + urand()*(borders[j].y - borders[j].x);
      evaluate(c);

      return c;
    }

    void initialize_population()
    {
      population.resize(population_size);
      best_population.resize(best_params_size);

      for (int i=0;i<population_size;i++)
        population[i] = initialize_creature();
    }

    void mutation(OptParams &params, float mutation_chance, float mutation_power)
    {
      Normal normal_gen = Normal(0, mutation_power);
      
      for (int i=0;i<MAX(1, mutation_chance * params.size());i++)
      {
        int id = urandi(0, params.size());
        float t = borders[id].y + 1;
        while (t >= borders[id].y || t <= borders[id].x)
          t = CLAMP(params[id] + normal_gen.get()*(borders[id].y - borders[id].x), borders[id].x, borders[id].y);
        params[id] = t;
      }
    }

    OptParams crossover(const OptParams &p1, const OptParams &p2)
    {
      int id = urandi(0, p1.size());
      OptParams p = p1;
      for (int i=id;i<p1.size();i++)
        p[i] = p2[i];
      return p;
    }

    void local_search(Creature &c, int iterations)
    {
      if (!c.local_optimizer)
      {
        UPGReconstructionResult res;
        res.structure = structure;
        res.parameters.p = c.params.differentiable;//TODO: fixme

        c.local_optimizer = get_optimizer_adam(func, local_opt_block, res);
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
          for (auto &v : c.params.differentiable)
            debug("%f ", v);
          debug("]\n");
        }
        best_population[worst_index] = c;
      }

      result_bins[CLAMP((int)(result_bin_count*c.loss), 0, result_bin_count-1)]++;
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
      std::sort(population.begin(), population.end(), [this](const Creature & a, const Creature & b) -> bool{    
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
    UPGFixedStructureMemeticOptimizer(UPGOptimizableFunction *_func, 
                                      const Block &settings, const UPGStructure &_structure) :
    UPGOptimizer(_func)
    {
      structure = _structure;
      gen = func->get_generator(structure);
      pd = func->get_full_parameters_description(gen.get());

      for (const auto &p : pd.get_block_params())
      {
        for (const auto &param_info : p.second.p)
        {
          if (param_info.type != ParameterType::CONST)
          {
            //logerr("added parameter %s", param_info.name.c_str());
            borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
          }
        }
      }

     verbose = settings.get_bool("verbose");
     finish_thr = settings.get_double("finish_threshold");
     budget = settings.get_int("iterations");
     local_opt_block.set_bool("verbose", settings.get_bool("verbose"));
     local_opt_block.set_bool("save_intermediate_images", false);
     local_opt_block.set_double("learning_rate", local_learning_rate);
    }
    virtual void optimize(int iters = -1) override
    {
      if (iters > 0)
        budget = iters;
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
          std::vector<float> chances = {0.1, 0.1 + 0.5, 0.1 + 0.5 + 0.4};
          float action_rnd = urand(0, chances.back());
          if (action_rnd < chances[0])
          {
            //create new random creature from scratch
            new_population[i] = initialize_creature();
          }
          else if (action_rnd < chances[1])
          {
            OptParams p = population[urandi(0, elites_count)].params;
            mutation(p, mutation_chance, urand()*mutation_power);
            new_population[i] = Creature(p, epoch);
            evaluate(new_population[i]);
          }
          else
          {
            auto parents = choose_parents_tournament(tournament_size);
            OptParams p = crossover(population[parents.first].params, population[parents.second].params);
            mutation(p, mutation_chance, urand()*mutation_power);
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
          debug("best res [");
          for (auto &c : best_population)
            debug("%.5f ", c.loss);
          debug("]\n");
        }

        int good_soulutions = 0;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr)
            good_soulutions++;

        float good_solutions_chance = (float)local_opt_count/good_soulutions;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr && (urand() < good_solutions_chance))
            local_search(population[i], local_opt_iters);

        if (verbose)
        {
          print_result_bins(result_bins);
          print_result_bins(current_bins);
        }

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
    int population_size = 2500;
    int best_params_size = 1;
    int budget = 200'000; //total number of function calls
    float mutation_chance = 0.1;
    float mutation_power = 0.1;
    int tournament_size = 128;
    int local_opt_count = 5;
    int local_opt_iters = 50;
    float local_learning_rate = 0.01;
    float good_soulution_thr = 0.05;
    float elites_fraction = 0.05;
    float finish_thr = 0;
    bool verbose = false;

    //generator-specific stuff
    UPGStructure structure;
    std::shared_ptr<UniversalGenInstance> gen;
    ParametersDescription pd;
    std::vector<glm::vec2> borders; //size equals total size of OptParams vector
    Block local_opt_block;

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
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_memetic(UPGOptimizableFunction *_func, 
                                                      const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGFixedStructureMemeticOptimizer>(_func, settings, _structure);
  }
}