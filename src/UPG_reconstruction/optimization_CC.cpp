#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>
namespace upg
{
  class UPGCCMemeticOptimizer : public UPGOptimizer
  {
  public:
    struct Creature;

    int population_size = 1000;
    int best_params_size = 1;
    int budget = 200'000; //total number of function calls
    int GA_generations = 5;
    int min_epochs = 10;
    float mutation_chance = 0.1;
    float mutation_power = 0.1;
    int tournament_size = 64;
    int local_opt_count = 2;
    int local_opt_iters = 50;
    float local_learning_rate = 0.01;
    float good_soulution_thr = 0.02;
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
    int parameters_count;

    int epoch = 0;
    int no_diff_function_calls = 0;
    int diff_function_calls = 0;

    //statistics
    constexpr static int result_bin_count = 1000;
    int result_bins[result_bin_count] = {};
    int current_bins[result_bin_count] = {};

    //creature is a unit of genetic algorithm
    //it is a set of parameters P with some additional info and 
    //values of loss(P) and fitness function
    struct Creature
    {
      Creature() = default;
      Creature(const OptParams &p)
      {
        params = p;
      }
      Creature(const OptParams &p, float l)
      {
        params = p;
        loss = l;
      }
      OptParams params;
      float loss = 1e9;

      OptParams best_params;
      float best_loss;
    };

    UPGCCMemeticOptimizer(UPGOptimizableFunction *_func, 
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
            logerr("added parameter %s", param_info.name.c_str());
            borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
          }
        }
      }
      parameters_count = borders.size();
      verbose = settings.get_bool("verbose");
      finish_thr = settings.get_double("finish_threshold");
      budget = settings.get_int("iterations");
      local_opt_block.set_bool("verbose", false);
      local_opt_block.set_bool("save_intermediate_images", false);
      local_opt_block.set_double("learning_rate", local_learning_rate);
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

    Creature initialize_creature()
    {
      Creature c;
      c.params.resize(borders.size());
      for (int j=0;j<borders.size();j++)
        c.params[j] = borders[j].x + urand()*(borders[j].y - borders[j].x);
      evaluate(c);
      c.best_params = c.params;
      c.best_loss = c.loss;

      return c;
    }

    void initialize_population()
    {
      population.resize(population_size);
      best_population.resize(best_params_size);

      for (int i=0;i<population_size;i++)
        population[i] = initialize_creature();
    }

    void mutation(OptParams &params, float mutation_chance, float mutation_power, std::pair<int, int> group)
    {
      Normal normal_gen = Normal(0, mutation_power);
      int group_size = group.second - group.first;
      
      for (int i=0;i<MAX(1, mutation_chance * group_size);i++)
      {
        int id = group.first + urandi(0, group_size);
        float t = borders[id].y + 1;
        while (t >= borders[id].y || t <= borders[id].x)
          t = CLAMP(params[id] + normal_gen.get()*(borders[id].y - borders[id].x), borders[id].x, borders[id].y);
        params[id] = t;
      }
    }

    OptParams crossover(const OptParams &p1, const OptParams &p2, std::pair<int, int> group)
    {
      int group_size = group.second - group.first;
      int id = group.first + urandi(0, group_size);
      OptParams p = p1;
      for (int i=id;i<p1.size();i++)
        p[i] = p2[i];
      //logerr("sizes %d %d %d", p1.size(), p2.size(), p.size());
      return p;
    }

    void local_search(Creature &c, int iterations)
    {
      UPGReconstructionResult res;
      res.structure = structure;
      res.parameters.p = func->opt_params_to_gen_params(c.params, pd);
      auto local_optimizer = get_optimizer_adam(func, local_opt_block, res);
      local_optimizer->optimize(iterations);
      UPGReconstructionResult lo_res = local_optimizer->get_best_results()[0];

      c.params = func->gen_params_to_opt_params(lo_res.parameters.p, pd);
      c.loss = lo_res.loss_optimizer;
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
    }

    void evaluate(Creature &c)
    {
      c.loss = func->f_no_grad(gen.get(), pd, c.params);
      register_new_param(c);
      no_diff_function_calls++;
      if (c.loss < c.best_loss)
      {
        c.best_loss = c.loss;
        c.best_params = c.params;
      }
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
      //for (int i=0;i<population_size;i++)
      // population[i].fitness = pos_fitness_function(i);
    }

    int parse_node_rec(const UPGStructure &structure, std::vector<std::pair<int, int>> &groups, bool merge_level, int start)
    {
      //how many children each node has. It should be done differently, but now it's just testing
      std::vector<int> node_powers       = {0, 0, 1, 2, 0, 0, 0, 0};
      std::vector<int> node_param_counts = {0, 1, 3, 0, 3, 2, 3, 1};
      int merge_node_num = 3;
      
      bool is_merge = merge_level && (structure.s[start] == merge_node_num);
      assert(structure.s[start] < node_powers.size());
      int power = node_powers[structure.s[start]];
      int pos = start+1;
      
      for (int i=0;i<power; i++)
      {
        int fin_pos = parse_node_rec(structure, groups, is_merge, pos);
        if (is_merge && structure.s[pos] != merge_node_num)
        {
          int p_st = groups.empty() ? 0 : groups.back().second;
          int p_cnt = 0;
          for (int p = pos; p < fin_pos; p++)
            p_cnt += node_param_counts[structure.s[p]];
          groups.push_back({p_st, p_st+p_cnt});
        }
        pos = fin_pos;
      }
      return pos;
    }

    std::vector<std::pair<int, int>> get_param_groups(const UPGStructure &structure)
    {
      std::vector<std::pair<int, int>> groups;
      parse_node_rec(structure, groups, true, 0);
      if (groups.empty())
        groups.push_back({0, borders.size()});
      for (auto &g : groups)
        logerr("group [%d %d]", g.first, g.second);
      return groups;
    }

    void perform_GA(int generations, std::pair<int, int> group)
    {
      for (int it = 0; it < generations; it++)
      {
        int elites_count = elites_fraction * population_size;
        std::vector<Creature> new_population(population_size);

        sort_and_calculate_fitness();
        for (int i = 0; i < elites_count; i++)
          new_population[i] = population[i];

        for (int i = elites_count; i < population_size; i++)
        {
          std::vector<float> chances = {0.1, 0.1 + 0.5, 0.1 + 0.5 + 0.4};
          float action_rnd = urand(0, chances.back());
          if (action_rnd < chances[0])
          {
            // create group parameters from scratch
            new_population[i].params = population[i].params;
            for (int p_n = group.first; p_n < group.second; p_n++)
              new_population[i].params[p_n] = borders[p_n].x + urand() * (borders[p_n].y - borders[p_n].x);
          }
          else if (action_rnd < chances[1])
          {
            OptParams p = population[urandi(0, elites_count)].params;
            mutation(p, mutation_chance, urand() * mutation_power, group);
            new_population[i] = Creature(p);
          }
          else
          {
            auto parents = choose_parents_tournament(tournament_size);
            OptParams p = crossover(population[parents.first].params, population[parents.second].params, group);
            mutation(p, mutation_chance, urand() * mutation_power, group);
            new_population[i] = Creature(p);
          }
          evaluate(new_population[i]);
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
        for (int i=0;i<local_opt_count;i++)
        {
          float rnd = urand(0, chances.back());
          for (int i=0;i<good_soulutions;i++)
          {
            if (rnd<chances[i])
            {
              local_search(population[i], local_opt_iters);
              break;
            }
          }
        }

        for (int i = 0; i < result_bin_count; i++)
          current_bins[i] = 0;
        for (auto &c : new_population)
          current_bins[CLAMP((int)(result_bin_count * c.loss), 0, result_bin_count - 1)]++;
        population = new_population;
      }
    }

    float dist(const Creature &c1, const Creature &c2)
    {
      float dist = 0;
      for (int i=0;i<borders.size();i++)
        dist += abs(c1.params[i]-c2.params[i])/abs(borders[i].x-borders[i].y);
      return dist/borders.size();
    }

    void perform_particle_swarm_optimization(int generations, std::pair<int, int> group)
    {
      float w = 0.729; // inertia weight
      float c1 = 1.49445; // cognitive weight
      float c2 = 1.49445; // social weight
      float r1, r2; // randomizers
      bool use_local_best = true;
      bool use_adaptive_weight = true;
      int group_size = group.second - group.first;

      std::vector<float> distances(population_size, 0);
      std::vector<int> indices(population_size, 0);
      std::vector<std::vector<float>> velocities(population_size, std::vector<float>(group_size,0));
      for (auto &v : velocities)
      {
        for (int j=0;j<group_size;j++)
          v[j] = 0.1*abs(borders[group.first+j].y - borders[group.first+j].x)*urand(-1,1);
      }

      for (int generation = 0; generation<generations; generation++)
      {        
        int success_count = 0;
        float total_movement = 0;
        Creature best_creature = best_population[0];

        for (int num=0;num<population_size;num++)
        {
          Creature &c = population[num];
          if (use_local_best)
          {
            for (int i=0;i<population_size;i++)
            {
              indices[i] = i;
              distances[i] = dist(c, population[i]);
            }

            std::sort(indices.begin(), indices.end(), [&distances](int a, int b)-> bool { return distances[a]<distances[b];});
            std::sort(indices.begin(), indices.begin() + population_size/5, [&](int a, int b)-> bool {return population[a].loss < population[b].loss;});

            best_creature = population[indices[0]];
          }

          OptParams best_params = best_creature.params;
          if (dist(c, best_creature) < 1e-5)
          {
            for (int j=0;j<velocities[num].size();j++)
              velocities[num][j] = 0;
            local_search(c, 1);
          }
          else
          {
            for (int j=0;j<velocities[num].size();j++)
            {
              r1 = urand();
              r2 = urand();
              velocities[num][j] = (w * velocities[num][j]) +
                              (c1 * r1 * (c.best_params[group.first+j] - c.params[group.first+j])) +
                              (c2 * r2 * (best_params[group.first+j] - c.params[group.first+j]));
            
              if (velocities[num][j] > abs(borders[group.first+j].x-borders[group.first+j].y))
                velocities[num][j] = abs(borders[group.first+j].x-borders[group.first+j].y);
              else if (velocities[num][j] < -abs(borders[group.first+j].x-borders[group.first+j].y))
                velocities[num][j] = -abs(borders[group.first+j].x-borders[group.first+j].y);
            
              c.params[group.first+j] = CLAMP(c.params[group.first+j]+velocities[num][j], borders[group.first+j].x, borders[group.first+j].y);
            }
            float prev_loss = c.best_loss;
            evaluate(c);
            if (c.loss < prev_loss)
              success_count++;
          }
          float vlen = 0;
          for (auto &v : velocities[num])
            vlen += SQR(v);
          total_movement += sqrt(vlen)/population_size;
        }

        //find worst particle and replace it with mutated version of the best one
        {
          int worst_pos = 0;
          for (int i=0;i<population_size;i++)
            if (population[i].best_loss > population[worst_pos].best_loss)
              worst_pos = i;
        
          population[worst_pos] = best_population[0];
          Normal normal_gen(0, borders.size()*(1-(no_diff_function_calls + diff_function_calls+0) / budget));
          int id = urandi(0, borders.size());
          float t = borders[id].y + 1;
          while (t >= borders[id].y || t <= borders[id].x)
            t = population[worst_pos].params[id] + normal_gen.get()*(borders[id].y - borders[id].x);
          population[worst_pos].params[id] = t;
        }

        float success_rate = (success_count+0.0)/population_size;
        if (use_adaptive_weight)
          w = success_rate;
        //logerr("TM SR %.3f %.3f", total_movement, success_rate);
        if (best_population[0].loss < finish_thr)
          break;
      }
    }

    virtual void optimize(int iters = -1) override
    {
      if (iters > 0)
        budget = iters;
      optimize_simple_CC();
    }

    void optimize_simple_CC()
    {      
      std::vector<std::pair<int, int>> groups = get_param_groups(structure);
      int predicted_epochs = budget/(groups.size()*GA_generations*(population_size+local_opt_count*local_opt_iters));
      while (predicted_epochs < min_epochs && groups.size() > 2)
      {
        debug("Expected number of epochs(%d) is too low. Cut groups count in half\n", predicted_epochs);
        std::vector<std::pair<int, int>> new_groups;
        for (int i=0;i<groups.size();i+=2)
        {
          if (i!=groups.size()-1)
            new_groups.push_back({groups[i].first, groups[i+1].second});
          else
            new_groups.push_back(groups[i]);
        }
        groups = new_groups;
        predicted_epochs = budget/(groups.size()*GA_generations*(population_size+local_opt_count*local_opt_iters));
      }
      std::vector<std::vector<std::vector<float>>> partial_solutions(groups.size());

      //Randomly initialize population, current best result is always best_population[0]
      if (population.empty())
        initialize_population();

      perform_GA(GA_generations, {0, borders.size()});
      
      //Get partial solutions to test later from our random population
      for (int g=0;g<groups.size();g++)
      {
        partial_solutions[g].resize(population.size());
        for (int i=0;i<population.size();i++)
        {
          partial_solutions[g][i].resize(groups[g].second - groups[g].first);
          for (int p_n=groups[g].first;p_n<groups[g].second;p_n++)
            partial_solutions[g][i][p_n-groups[g].first] = population[i].params[p_n];
        }
      }

      while (no_diff_function_calls + diff_function_calls < budget)
      {
        //Optimize values group by group
        for (int g=0;g<groups.size();g++)
        {
          //Set stored partial solutions to our group and best parameters for other groups
          for (int i=0;i<population_size;i++)
          {
            population[i] = best_population[0];
            for (int p_n=groups[g].first;p_n<groups[g].second;p_n++)
              population[i].params[p_n] = partial_solutions[g][i][p_n-groups[g].first];
            evaluate(population[i]);
          }

          //Perform standart genetic algorithm
          perform_GA(GA_generations, groups[g]);

          //save improved parameters for this group to partial solution storage
          for (int i=0;i<population.size();i++)
          {
            for (int p_n=groups[g].first;p_n<groups[g].second;p_n++)
              partial_solutions[g][i][p_n-groups[g].first] = population[i].params[p_n];
          }
        }

        //local optimization for the whole parameters' set
        //Creature c = best_population[0];
        //local_search(c, local_opt_iters);


        if (verbose)
        {
          debug("EPOCH %d (%d+%d/%d)\n", epoch, no_diff_function_calls, diff_function_calls, budget);
          debug("best res [");
          for (auto &c : best_population)
            debug("%.6f ", c.loss);
          debug("]\n");
        }


        /*
        int good_soulutions = 0;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr)
            good_soulutions++;

        float good_solutions_chance = (float)local_opt_count/good_soulutions;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr && (urand() < good_solutions_chance))
            local_search(population[i], local_opt_iters);
        */

        if (verbose)
        {
          //print_result_bins(result_bins);
          //print_result_bins(current_bins);
        }

        epoch++;
        if (best_population[0].loss < finish_thr)
          break;
      }
    }
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_CC(UPGOptimizableFunction *_func,
                                                 const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGCCMemeticOptimizer>(_func, settings,_structure);
  }
}