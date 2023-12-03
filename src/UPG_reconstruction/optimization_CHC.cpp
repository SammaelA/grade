#include "optimization.h"
#include "common_utils/distribution.h"
#include <algorithm>
namespace upg
{
  class UPGCHCMemeticOptimizer : public UPGOptimizer
  {
  public:
    struct Creature;

    int population_size = 250;
    int best_params_size = 1;
    int budget = 150'000; //total number of function calls
    int local_opt_count = 5;
    int local_opt_iters = 50;
    float local_learning_rate = 0.01;
    float good_soulution_thr = 0.05;
    float divergence_rate = 0.3;
    float distance_thr = 0.01;

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
      int birth_epoch = -1;
      int local_opt_num = 0;
    };

    UPGCHCMemeticOptimizer(UPGOptimizableFunction *_func, 
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
      local_opt_block.set_bool("verbose", settings.get_bool("verbose"));
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
        debug("%d new best (%f) [", no_diff_function_calls, c.loss);
        for (auto &v : c.params.differentiable)
          debug("%f ", v);
        debug("]\n");
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

    float initialize_CHC_distance(bool initial_dist)
    {
      if (initial_dist)
        return population_size/4.0;
      else 
        return divergence_rate*(1-divergence_rate)*population_size;
    }

    std::vector<std::pair<int,int>> get_pairs(std::vector<Creature> population)
    {
      std::vector<int> indices(population.size(), 0);
      for (int i=0;i<population.size(); i++)
        indices[i] = i;
      std::random_shuffle(indices.begin(),indices.end());
      std::vector<std::pair<int,int>> pairs(population.size()/2);
      for (int i=0;i<population.size(); i+=2)
        pairs[i/2] = {indices[i], indices[i+1]};
      return pairs;
    }

    float distance(const Creature &c1, const Creature &c2)
    {
      float dist = 0;
      for (int i=0;i<parameters_count;i++)
        dist += abs(c1.params[i] - c2.params[i])/(borders[i].y - borders[i].x);
      return dist;
    }

    void recombine(Creature &c1, Creature &c2)
    {
      //swap half the differing params at random
      for (int i=0;i<parameters_count;i++)
      {
        if (abs(c1.params[i] - c2.params[i])/(borders[i].y - borders[i].x) > distance_thr && urand() > 0.5)
          std::swap(c1.params[i], c2.params[i]);
      }
      evaluate(c1);
      evaluate(c2);
    }

    std::vector<Creature> recombine(const std::vector<Creature> &prev_pop, float CHC_distance)
    {
      std::vector<Creature> add_pop;
      add_pop.reserve(prev_pop.size());
      auto pairs = get_pairs(prev_pop);

      for (auto p : pairs)
      {
        //logerr("distance %f", distance(prev_pop[p.first], prev_pop[p.second]));
        if (0.25*distance(prev_pop[p.first], prev_pop[p.second]) > CHC_distance)
        {
          add_pop.push_back(prev_pop[p.first]);
          add_pop.push_back(prev_pop[p.second]);
          recombine(add_pop[add_pop.size()-2], add_pop[add_pop.size()-1]);
        }
      }
      logerr("%d/%d pairs recombined", (int)(add_pop.size()/2), (int)(pairs.size()));

      return add_pop;
    }

    std::vector<Creature> select(const std::vector<Creature> &prev_pop, std::vector<Creature> &add_pop, bool *equals)
    {
      std::vector<Creature> new_pop = prev_pop;
      auto cmp_f = [](const Creature & a, const Creature & b) -> bool {return a.loss < b.loss;};
      std::sort(new_pop.begin(), new_pop.end(), cmp_f);
      std::sort(add_pop.begin(), add_pop.end(), cmp_f);
      int i = 0;
      int j = new_pop.size()-1;
      while (i<add_pop.size() && j >= 0 && new_pop[j].loss >= add_pop[i].loss)
      {
        new_pop[j] = add_pop[i];
        i++;
        j--;
      }
      *equals = (i==0);
      return new_pop;
    }

    void diverge(std::vector<Creature> &pop, const std::vector<Creature> &prev_pop)
    {
      auto cmp_f = [](const Creature & a, const Creature & b) -> bool {return a.loss < b.loss;};
      Creature c = *(std::max_element(prev_pop.begin(), prev_pop.end(), cmp_f));
      pop[0] = c;

      for (int i=0;i<pop.size();i++)
      {
        pop[i] = c;
        for (int j=0;j<c.params.size();j++)
          if (urand() < divergence_rate)
            pop[i].params[j] = borders[j].x + urand()*(borders[j].y - borders[j].x);
        evaluate(pop[i]);
      }
    }
    
    virtual void optimize(int iters = -1) override
    {
      initialize_population();
      float CHC_distance = initialize_CHC_distance(true);
      int it = 0;
      while (no_diff_function_calls + diff_function_calls < budget)
      {
        std::vector<Creature> add_pop = recombine(population, CHC_distance);
        bool equals;
        std::vector<Creature> new_pop = select(population, add_pop, &equals);
        if (equals)
          CHC_distance--;
        if (CHC_distance <= 0)
        {
          logerr("diverge");
          diverge(new_pop, population);
          CHC_distance = initialize_CHC_distance(false);
        }
        population = new_pop;
        if (it % 25 == 0)
          print_result_bins(result_bins);
        it++;
      }
    }
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_CHC(UPGOptimizableFunction *_func,
                                                  const Block &settings, const UPGStructure &_structure)
  {
    return std::make_shared<UPGCHCMemeticOptimizer>(_func, settings,_structure);
  }
}