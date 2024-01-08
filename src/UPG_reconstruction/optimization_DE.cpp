#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>

namespace upg
{
  //Adaptive weight taken from paper "A novel particle swarm optimization algorithm with adaptive inertia weight"
  class UPGDifferentiableEvolutionOptimizer : public UPGOptimizer
  {
  private:
    //creature is a unit of genetic algorithm
    //it is a set of parameters P with some additional info and 
    //values of loss(P) and fitness function
    struct Creature
    {
      OptParams params;
      float loss;
      float F = 0.5;
      float CR = 0.5;
      
      std::shared_ptr<UPGOptimizer> local_optimizer;
      int local_opt_num = 0;
    };

    Creature initialize_creature()
    {
      Creature c;
      c.params.resize(borders.size());
      for (int j=0;j<borders.size();j++)
        c.params[j] = borders[j].x + urand()*(borders[j].y - borders[j].x);
      c.F = urand();
      c.CR = urand();
      evaluate(c);

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
          for (int i=0;i<c.params.size();i++)
            debug("%f ", c.params[i]);
          debug("]\n");
        }
        best_population[worst_index] = c;
      }

      result_bins[CLAMP((int)(result_bin_count*c.loss), 0, result_bin_count-1)]++;
    
      double slide = 0.9975;
      for (int i=0;i<borders.size();i++)
      {
        int bin = stat[i].size()*((c.params[i] - borders[i].x)/(borders[i].y-borders[i].x));
        bin = CLAMP(bin, 0, stat[i].size()-1);
        stat[i][bin].first++;
        stat[i][bin].second = slide*stat[i][bin].second + (1-slide)*c.loss;
      }
    }

    void evaluate(Creature &c)
    {
      c.loss = func->f_no_grad(gen.get(), pd, c.params);
      register_new_param(c);
      no_diff_function_calls++;
    }

    void local_search(Creature &c, int iterations)
    {
      if (!c.local_optimizer)
      {
        UPGReconstructionResult res;
        res.structure = structure;
        res.parameters.p = func->opt_params_to_gen_params(c.params, pd);

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
    UPGDifferentiableEvolutionOptimizer(UPGOptimizableFunction *_func, 
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
      stat.resize(borders.size());
      for (auto &a : stat)
        for (auto &p : a)
          p = {0, 0.0};

      verbose = settings.get_bool("verbose");
      finish_thr = settings.get_double("finish_threshold");
      budget = settings.get_int("iterations");
      local_opt_block.set_bool("verbose", false);
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

    OptParams rand_1_mutation(int i, float F)
    {
      OptParams v = population[i].params;
      int r1=i,r2=i,r3=i;
      while (r1==i)
        r1 = urandi(0, population_size);
      while (r2==i || r2==r1)
        r2 = urandi(0, population_size);
      while (r3==i || r3==r1 || r3==r2)
        r3 = urandi(0, population_size);
      
      for (int j=0;j<v.size();j++)
        v[j] = CLAMP(best_population[0].params[j] + F*(population[r2].params[j] - population[r3].params[j]), borders[j].x, borders[j].y);
      
      return v;
    }

    OptParams binary_crossover(const OptParams &x, const OptParams &v, float CR)
    {
      OptParams u = x;
      int j_rand = urandi(0, v.size());
      for (int j=0;j<v.size();j++)
      {
        if (urand() < CR)
          u[j] = v[j];
      }
      u[j_rand] = v[j_rand];

      return u;
    }

    float get_new_F(float F)
    {
      float Fl = 0.1, Fu = 0.9, tau_1 = 0.1;
      if (urand() < tau_1)
        return Fl + urand()*Fu;
      else
        return F;
    }

    float get_new_CR(float CR)
    {
      float tau_2 = 0.1;
      if (urand() < tau_2)
        return urand();
      else
        return CR;
    }

    virtual void optimize(int iters = -1) override
    {
      population_size = 250;
      int l = 1;
      int k = 50;
      float expected_replace_rate = 0.25;

      if (iters > 0)
        budget = iters;
      if (population.empty())
        initialize_population();

      while (no_diff_function_calls + diff_function_calls < budget)
      {        
        int replaced = 0;
        for (int i=0;i<population_size;i++)
        {
          float new_F = get_new_F(population[i].F);
          float new_CR = get_new_CR(population[i].CR);
          OptParams x = population[i].params;
          OptParams v = rand_1_mutation(i, new_F);
          OptParams u = binary_crossover(x, v, new_CR);
          Creature new_c;
          new_c.params = u;
          new_c.F = new_F;
          new_c.CR = new_CR;
          evaluate(new_c);
          if (new_c.loss < population[i].loss)
          {
            replaced++;
            population[i] = new_c;
          }
        }
        if (epoch > 0 && epoch % l == 0)
        {
          std::sort(population.begin(), population.end(), [](const Creature & a, const Creature & b) -> bool {    
                    return a.loss < b.loss;});
          for (int i=0;i<population_size-k;i++)
          {
            if (urand() < 1/(i+1.0))
              local_search(population[i], 5);
          }
          for (int i=population_size-k;i<population_size;i++)
            population[i] = initialize_creature();
        }
        if (verbose)
        {
          debug("EPOCH %d (%d %d): replaced %d/%d (k=%d). Best %.8f\n", epoch, no_diff_function_calls, diff_function_calls,
                                                                        replaced, population_size, k, best_population[0].loss);
          for (auto &arr : stat)
          {
            //debug("param stat: ");
            //for (auto &p : arr)
            //  debug("[%f %d]",(float)p.second, p.first);
            //debug("\n");
          }
        }
        
        if (replaced < expected_replace_rate*population_size)
          k++;
        else
          k--;
        k = CLAMP(k, 0.05*population_size, 0.75*population_size);

        if (best_population[0].loss < finish_thr)
          break;
        epoch++;
      }
      if (verbose)
      {
        for (auto &p : population)
          logerr("[%f %f] %f",p.F, p.CR, p.loss);
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
    int population_size = 30;
    int best_params_size = 1;
    int budget = 200'000; //total number of function calls
    float local_learning_rate = 0.01;
    float finish_thr = 0;
    bool verbose = false;

    //generator-specific stuff
    UPGStructure structure;
    std::shared_ptr<UniversalGenInstance> gen;
    ParametersDescription pd;
    std::vector<glm::vec2> borders; //size equals total size of OptParams vector
    std::vector<std::array<std::pair<int, double>, 8>> stat;
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
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_differentiable_evolution(UPGOptimizableFunction *_func, 
                                                                       const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGDifferentiableEvolutionOptimizer>(_func, settings, _structure);
  }
}