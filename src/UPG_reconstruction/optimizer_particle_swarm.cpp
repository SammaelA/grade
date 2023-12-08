#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>

namespace upg
{
  class UPGParticleSwarmOptimizer : public UPGOptimizer
  {
  private:
    //creature is a unit of genetic algorithm
    //it is a set of parameters P with some additional info and 
    //values of loss(P) and fitness function
    struct Creature
    {
      OptParams params;
      float loss;
      std::vector<float> velocity;

      OptParams best_params;
      float best_loss;
      
      std::shared_ptr<UPGOptimizer> local_optimizer;
      int local_opt_num = 0;
    };

    Creature initialize_creature()
    {
      Creature c;
      c.params.resize(borders.size());
      c.velocity.resize(borders.size());

      for (int j=0;j<borders.size();j++)
      {
        c.params[j] = borders[j].x + urand()*(borders[j].y - borders[j].x);
        c.velocity[j] = 0.1*abs(borders[j].y - borders[j].x)*urand(-1,1);
      }
      evaluate(c);

      c.best_params = c.params;
      c.best_loss = c.loss;

      return c;
    }

    void initialize_population()
    {
      population.resize(population_size);
      best_population.resize(best_params_size);
      for (auto &bp : best_population)
        bp.loss = 1e9;

      for (int i=0;i<population_size;i++)
        population[i] = initialize_creature();
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
      if (c.loss < c.best_loss)
      {
        c.best_loss = c.loss;
        c.best_params = c.params;
      }
      register_new_param(c);
      no_diff_function_calls++;
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
    UPGParticleSwarmOptimizer(UPGOptimizableFunction *_func, 
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
     local_opt_block.set_bool("verbose", settings.get_bool("verbose"));
     local_opt_block.set_bool("save_intermediate_images", false);
     local_opt_block.set_double("learning_rate", local_learning_rate);
    }

    float dist(const Creature &c1, const Creature &c2)
    {
      float dist = 0;
      for (int i=0;i<borders.size();i++)
        dist += abs(c1.params[i]-c2.params[i])/abs(borders[i].x-borders[i].y);
      return dist/borders.size();
    }

    virtual void optimize(int iters = -1) override
    {
      if (iters > 0)
        budget = iters;
      if (population.empty())
        initialize_population();

      float w = 0.729; // inertia weight
      float c1 = 1.49445; // cognitive weight
      float c2 = 1.49445; // social weight
      float r1, r2; // randomizers

      std::vector<float> distances(population_size, 0);
      std::vector<int> indices(population_size, 0);

      while (no_diff_function_calls + diff_function_calls < budget)
      {        
        float total_movement = 0;
        for (auto &c : population)
        {
          for (int i=0;i<population_size;i++)
          {
            indices[i] = i;
            distances[i] = dist(c, population[i]);
          }

          std::sort(indices.begin(), indices.end(), [&distances](int a, int b)-> bool { return distances[a]<distances[b];});
          std::sort(indices.begin(), indices.begin() + population_size/5, [&](int a, int b)-> bool {return population[a].loss < population[b].loss;});

          OptParams local_best = population[indices[0]].params;
          if (dist(c, population[indices[0]]) < 1e-5)
          {
            for (int j=0;j<c.velocity.size();j++)
              c.velocity[j] = 0;
            local_search(c, 1);
          }
          else
          {
            for (int j=0;j<c.velocity.size();j++)
            {
              r1 = urand();
              r2 = urand();
              c.velocity[j] = (w * c.velocity[j]) +
                              (c1 * r1 * (c.best_params[j] - c.params[j])) +
                              (c2 * r2 * (local_best[j] - c.params[j]));
            
              if (c.velocity[j] > abs(borders[j].x-borders[j].y))
                c.velocity[j] = abs(borders[j].x-borders[j].y);
              else if (c.velocity[j] < -abs(borders[j].x-borders[j].y))
                c.velocity[j] = -abs(borders[j].x-borders[j].y);
            
              c.params[j] = CLAMP(c.params[j]+c.velocity[j], borders[j].x, borders[j].y);
            }
            evaluate(c);
          }
          float vlen = 0;
          for (auto &v : c.velocity)
            vlen += SQR(v);
          total_movement += sqrt(vlen)/population_size;
        }
        if (verbose)
        {
          //print_result_bins(result_bins);
          //print_result_bins(current_bins);
        }
        logerr("total movement %f", total_movement);

        epoch++;
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
    int population_size = 1000;
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

  std::shared_ptr<UPGOptimizer> get_optimizer_particle_swarm(UPGOptimizableFunction *_func, 
                                                             const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGParticleSwarmOptimizer>(_func, settings, _structure);
  }
}