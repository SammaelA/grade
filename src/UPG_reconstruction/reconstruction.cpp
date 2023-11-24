#include "upg.h"
#include "preprocessing.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"
#include "tinyEngine/engine.h"
#include "custom_diff_render/custom_diff_render.h"
#include "simple_render_and_compare.h"
#include "graphics_utils/modeling.h"
#include "common_utils/distribution.h"
#include <memory>
#include <unistd.h>
#include <algorithm>

namespace upg
{
  std::vector<ParametersDescription::Param> get_camera_params(const CameraSettings &cam)
  {
    std::vector<ParametersDescription::Param> params;

    params.push_back({cam.origin.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.x"});
    params.push_back({cam.origin.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.y"});
    params.push_back({cam.origin.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.z"});

    params.push_back({cam.target.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.x"});
    params.push_back({cam.target.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.y"});
    params.push_back({cam.target.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.z"});

    params.push_back({cam.up.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.x"});
    params.push_back({cam.up.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.y"});
    params.push_back({cam.up.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.z"});

    params.push_back({cam.z_near, 0.001f, 1.0f, ParameterType::MUTABLE_FLOAT, "z_near"});
    params.push_back({cam.z_far, 10.0f, 1000.0f, ParameterType::MUTABLE_FLOAT, "z_far"});
    params.push_back({cam.fov_rad, 0.1f, 1.5f, ParameterType::MUTABLE_FLOAT, "fov_rad"});

    return params;
  }

  unsigned get_camera_block_id(int camera_n)
  {
    return (1 << 16) + (unsigned)camera_n;
  }

  ParametersDescription get_cameras_parameter_description(const std::vector<ReferenceView> &reference_views)
  {
    ParametersDescription pd;
    for (int i = 0; i < reference_views.size(); i++)
    {
      auto pv = get_camera_params(reference_views[i].camera);
      if (reference_views[i].fixed_camera)
      {
        for (auto &p : pv)
          p.type = ParameterType::CONST;
      }
      pd.add_parameters(get_camera_block_id(i), "camera_" + std::to_string(i), pv);
    }
    return pd;
  }

  CameraSettings camera_from_params(const std::vector<float> &p)
  {
    CameraSettings cam;
    cam.origin = glm::vec3(p[0], p[1], p[2]);
    cam.target = glm::vec3(p[3], p[4], p[5]);
    cam.up     = glm::vec3(p[6], p[7], p[8]);
    cam.z_near = p[9];
    cam.z_far = p[10];
    cam.fov_rad = p[11];
    return cam;
  }

  
  void grad_jac_mult(const UniversalGenJacobian &dPos_dP, std::span<const float> dLoss_dPos, 
                     std::span<float> out_dLoss_dP)
  {
    assert(out_dLoss_dP.size() == dPos_dP.get_yn());
    assert(dLoss_dPos.size() == dPos_dP.get_xn());
    for (int i=0;i<dPos_dP.get_yn();i++)
    { out_dLoss_dP[i] = 0;
      for (int j=0;j<dPos_dP.get_xn();j++)
        out_dLoss_dP[i] += dPos_dP.at(i,j)*dLoss_dPos[j]; 
    }
  }

  class UPGOptimizer
  {
  public:
    virtual ~UPGOptimizer() = default;
    virtual void optimize(int iters = -1) = 0; //iters = -1 mean that we take it from settings
    virtual std::vector<UPGReconstructionResult> get_best_results() = 0;
    UPGOptimizer(IDiffRender *_diff_render, IDiffRender *_simple_render, const ParametersDescription &_cameras_pd)
    {
      diff_render = _diff_render;
      simple_render = _simple_render;
      cameras_pd = _cameras_pd;
    }
  protected:
    //all parameters that can be changed by optimizer structured
    //in a convenient (for optimizer) way
    struct Params
    {
      std::vector<float> differentiable;
      int size() const
      {
        return differentiable.size();
      }
      void resize(int size)
      {
        differentiable.resize(size);
      }
      float &operator[](int index)
      {
        return differentiable[index];
      }
      const float &operator[](int index) const
      {
        return differentiable[index];
      }
    };

    void opt_params_to_gen_params_and_camera(const Params &params, const ParametersDescription &pd, 
                                            /*out*/ std::vector<float> &full_gen_params,
                                            /*out*/ std::vector<CameraSettings> &cameras)
    {
      //iterate all parameters' groups from description
      //map orders them by block_id, so generator's params are first, and cameras' after it
      int diff_i = 0;
      full_gen_params.reserve(pd.get_total_params_count());
      for (const auto &p : pd.get_block_params())
      {
        std::vector<float> camera_params;
        bool is_camera_block = p.first >= get_camera_block_id(0);
        std::vector<float> &p_v = is_camera_block ? camera_params : full_gen_params;
        for (auto &par : p.second.p)
        {
          if (par.type == ParameterType::CONST)
            p_v.push_back(par.value);
          else if (par.type == ParameterType::DIFFERENTIABLE)
          {
            float v = CLAMP(params.differentiable[diff_i], par.min_val, par.max_val);
            p_v.push_back(v);
            diff_i++;
          }
          else
          {
            // TODO: other types
          }
        }
        if (is_camera_block)
          cameras.push_back(camera_from_params(camera_params));
      }
    }

    //calculate function that we optimize and it's gradient (put into given span)
    //requires already created UniversalGenInstance
    //size of out_grad - is a number of differentiable parameters in params
    float f_grad_f(UniversalGenInstance &gen, const ParametersDescription &pd,
                   const Params &params, std::span<float> out_grad)
    {
      std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
      std::vector<CameraSettings> cameras;
      opt_params_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenJacobian dPos_dP;
      UniversalGenMesh mesh = gen.generate(full_gen_params, &dPos_dP);
      float res = diff_render->render_and_compare_silhouette(mesh.pos, cameras);
      std::span<const float> dLoss_dPos(diff_render->get_vertex_grad(), mesh.pos.size());
      grad_jac_mult(dPos_dP, dLoss_dPos, out_grad);

      return res;
    }

    float f_no_grad(UniversalGenInstance &gen, const ParametersDescription &pd, const Params &params)
    {
      std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
      std::vector<CameraSettings> cameras;
      opt_params_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenMesh mesh = gen.generate(full_gen_params, nullptr);
      
      return simple_render->render_and_compare_silhouette(mesh.pos, cameras);
    }

    int render_w, render_h;
    IDiffRender *diff_render;
    IDiffRender *simple_render;
    ParametersDescription cameras_pd;
  };

  class UPGOptimizerAdam : public UPGOptimizer
  {
  public:
    UPGOptimizerAdam(IDiffRender *_diff_render, IDiffRender *_simple_render, const ParametersDescription &_cameras_pd, 
                     const Block &settings, const UPGReconstructionResult &start_params) :
    UPGOptimizer(_diff_render, _simple_render, _cameras_pd),
    gen(start_params.structure)
    {
      iterations = settings.get_int("iterations");
      alpha = settings.get_double("learning_rate", 0.01);
      beta_1 = settings.get_double("beta_1", 0.9);
      beta_2 = settings.get_double("beta_2", 0.999);
      eps = settings.get_double("eps", 1e-8);
      verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;

      X_n = start_params.parameters.p.size();
      X.differentiable = start_params.parameters.p;//TODO: fixme
      gen_structure = start_params.structure;

      pd.add(cameras_pd);
      pd.add(gen.desc);

      V = std::vector<float>(X_n, 0); 
      S = std::vector<float>(X_n, 0);
      best_params = X;
      x_grad = std::vector<float>(X_n, 0); 
      best_result = 1e9;
    }
    virtual void optimize(int iters = -1) override
    {
      for (int iter=0; iter< (iters>0 ? iters : iterations); iter++)
      {
        float val = f_grad_f(gen, pd, X, x_grad);
        if (val < best_result)
        {
          best_params = X;
          best_result = val;
        }
        for (int i=0;i<X_n;i++)
        {
          float g = x_grad[i];
          V[i] = beta_1 * V[i] + (1-beta_1)*g;
          float Vh = V[i] / (1 - pow(beta_1, total_iterations+1)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
          float Sh = S[i] / (1 - pow(beta_2, total_iterations+1)); 
          X.differentiable[i] -= alpha*Vh/(sqrt(Sh) + eps);
        }
        if ((total_iterations % 5 == 0) && verbose)
          debug("Adam iter %3d  val = %.6f best_val = %.6f\n", total_iterations, val, best_result);
        total_iterations++;
      }
      if (verbose)
        debug("Adam final res val = %.6f best_val = %.6f\n", best_result, best_result);
    }

    virtual std::vector<UPGReconstructionResult> get_best_results() override
    {
      UPGReconstructionResult res;
      res.structure = gen_structure;
      {
        std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
        std::vector<CameraSettings> cameras;
        opt_params_to_gen_params_and_camera(best_params, pd, full_gen_params, cameras);
        res.parameters.p = full_gen_params;
        res.loss_optimizer = best_result;
      }
      return {res};
    }
    Params get_best_params_in_optimizer_format() { return best_params;}
    float get_best_result_in_optimizer_format() { return best_result;}

  private:
    int iterations;
    float alpha, beta_1, beta_2, eps;
    bool verbose;
    int X_n;
    UPGOptimizer::Params X;
    UniversalGenInstance gen;
    ParametersDescription pd;
    UPGStructure gen_structure;

    std::vector<float> V; 
    std::vector<float> S;
    UPGOptimizer::Params best_params;
    std::vector<float> x_grad; 
    float best_result;
    int total_iterations = 0;
  };

  class UPGFixedStructureMemeticOptimizer : public UPGOptimizer
  {
  private:
    //creature is a unit of genetic algorithm
    //it is a set of parameters P with some additional info and 
    //values of loss(P) and fitness function
    struct Creature
    {
      Creature() = default;
      Creature(const Params &p, int epoch)
      {
        params = p;
        birth_epoch = epoch;
        local_opt_num = 0;
      }
      Creature(const Params &p, float l, int epoch)
      {
        params = p;
        loss = l;
        birth_epoch = epoch;
        local_opt_num = 0;
      }
      Params params;
      std::shared_ptr<UPGOptimizerAdam> local_optimizer;
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

    void mutation(Params &params, float mutation_chance, float mutation_power)
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

    Params crossover(const Params &p1, const Params &p2)
    {
      int id = urandi(0, p1.size());
      Params p = p1;
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

        c.local_optimizer = std::make_shared<UPGOptimizerAdam>(diff_render, simple_render, cameras_pd, local_opt_block, res);
      }
      c.local_optimizer->optimize(iterations);

      c.params = c.local_optimizer->get_best_params_in_optimizer_format();
      c.loss = c.local_optimizer->get_best_result_in_optimizer_format();
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
        for (auto &v : c.params.differentiable)//TODO: fixme
          debug("%f ", v);
        debug("]\n");
        best_population[worst_index] = c;
      }

      result_bins[CLAMP((int)(result_bin_count*c.loss), 0, result_bin_count-1)]++;
    }

    void evaluate(Creature &c)
    {
      c.loss = f_no_grad(gen, pd, c.params);
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
    UPGFixedStructureMemeticOptimizer(IDiffRender *_diff_render, IDiffRender *_simple_render, const ParametersDescription &_cameras_pd, 
                                      const Block &settings, const UPGStructure &_structure) :
    UPGOptimizer(_diff_render, _simple_render, _cameras_pd),
    gen(_structure)
    {
      structure = _structure;
      pd.add(cameras_pd);
      pd.add(gen.desc);

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

     local_opt_block.set_int("render_w", render_w);
     local_opt_block.set_int("render_h", render_h);
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
            Params p = population[urandi(0, elites_count)].params;
            mutation(p, mutation_chance, urand()*mutation_power);
            new_population[i] = Creature(p, epoch);
            evaluate(new_population[i]);
          }
          else
          {
            auto parents = choose_parents_tournament(tournament_size);
            Params p = crossover(population[parents.first].params, population[parents.second].params);
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

        debug("EPOCH %d (%d+%d/%d)\n", epoch, no_diff_function_calls, diff_function_calls, budget);
        debug("best res [");
        for (auto &c : best_population)
          debug("%.5f ", c.loss);
        debug("]\n");

        int good_soulutions = 0;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr)
            good_soulutions++;

        float good_solutions_chance = (float)local_opt_count/good_soulutions;
        for (int i=0;i<population_size;i++)
          if (population[i].loss < good_soulution_thr && (urand() < good_solutions_chance))
            local_search(population[i], local_opt_iters);

        print_result_bins(result_bins);
        print_result_bins(current_bins);
        //debug("BINS\n");
        //for (int i=0;i<100;i++)
        //  debug("%d) %d\n",i, result_bins[i]);
        epoch++;
      }
    }
    virtual std::vector<UPGReconstructionResult> get_best_results() override
    {
      UPGReconstructionResult res;
      res.structure = structure;
      {
        std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
        std::vector<CameraSettings> cameras;
        opt_params_to_gen_params_and_camera(best_population[0].params, pd, full_gen_params, cameras);
        res.parameters.p = full_gen_params;
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
    float good_soulution_thr = 0.015;
    float elites_fraction = 0.05;

    //generator-specific stuff
    UPGStructure structure;
    UniversalGenInstance gen;
    ParametersDescription pd;
    std::vector<glm::vec2> borders; //size equals total size of Params vector
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

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk)
  {
    //load settings from given blk
    Block *input_blk = blk.get_block("input");
    Block *gen_blk   = blk.get_block("generator");
    Block *opt_blk   = blk.get_block("optimization");
    Block *res_blk   = blk.get_block("results");
    if (!input_blk || !gen_blk || !opt_blk || !res_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    //get ReconstructionReference - all info about the object that we want to reconstruct
    ReconstructionReference reference = get_reference(*input_blk);

    //get start parameters for optimization. They are required for Adam and other local optimizers
    //and have to be set manually
    UPGReconstructionResult start_params;
    Block *start_params_blk = opt_blk->get_block("start");
    if (start_params_blk)
    {
      start_params_blk->get_arr("params", start_params.parameters.p);
      start_params_blk->get_arr("structure", start_params.structure.s);
    }

    //perform optimization. There might be one or several steps of it, I expect the first step 
    //to be some sort of Genetic Algorithm and others - Adam optimizers for fine-tuning the params
    int step_n = 0;
    std::vector<UPGReconstructionResult> opt_res = {start_params};
    while (opt_blk->get_block("step_"+std::to_string(step_n)))
    {
      Block *step_blk = opt_blk->get_block("step_"+std::to_string(step_n));

      int render_w = step_blk->get_int("render_w", 128);
      int render_h = step_blk->get_int("render_h", 128);
      auto cameras_pd = get_cameras_parameter_description(reference.images);
      std::unique_ptr<IDiffRender> diff_render(get_halfgpu_custom_diff_render());
      std::unique_ptr<IDiffRender> simple_render(new NonDiffRender());
      
      //get diff_render settings
      IDiffRender::Settings diff_render_settings;
      diff_render_settings.image_w = render_w;
      diff_render_settings.image_h = render_h;

      std::vector<Texture> references;
      for (int i=0; i<reference.images.size(); i++)
      {
        reference.images[i].resized_mask = resize_mask(reference.images[i].mask, render_w, render_h, true);
        references.push_back(reference.images[i].resized_mask);
      }

      if (diff_render)
      diff_render->init_optimization(references, diff_render_settings, 
                                     step_blk->get_bool("save_intermediate_images", false));
      if (simple_render)
      simple_render->init_optimization(references, diff_render_settings, 
                                       step_blk->get_bool("save_intermediate_images", false));

      std::unique_ptr<UPGOptimizer> optimizer;
      std::string optimizer_name = step_blk->get_string("optimizer_name", "adam");
      if (optimizer_name == "adam")
        optimizer.reset(new UPGOptimizerAdam(diff_render.get(), simple_render.get(), cameras_pd, *step_blk, opt_res[0]));
      else if (optimizer_name == "memetic")
        optimizer.reset(new UPGFixedStructureMemeticOptimizer(diff_render.get(), simple_render.get(), cameras_pd, *step_blk, start_params.structure));
      optimizer->optimize();
      opt_res = optimizer->get_best_results();
      step_n++;
    }
    
    //TODO: compare results with reference and calculate reconstruction quality
    //Also calculate another quality metric if we have synthetic reference
    for (auto &result : opt_res)
    {
      ComplexModel reconstructed_model;
      std::string save_directory = "saves/" + res_blk->get_string("save_folder") + "/";
      if (res_blk->get_string("save_folder") != "")
        prepare_directory("saves/" + res_blk->get_string("save_folder"));
      if (!create_model(result.structure, result.parameters, reconstructed_model))
        logerr("failed to create model from reconstructed structure and parameters!");
      if (res_blk->get_bool("check_image_quality"))
        result.quality_ir = get_image_based_quality(reference, reconstructed_model);
      if (reference.model.is_valid() && res_blk->get_bool("check_model_quality"))
        result.quality_synt = get_model_based_quality(reference, reconstructed_model);
    
      if (res_blk->get_bool("save_turntable") && res_blk->get_block("save_turntable_hydra_settings"))
        render_model_turntable(*(res_blk->get_block("save_turntable_hydra_settings")), reconstructed_model);
      if (reference.model.is_valid() && res_blk->get_bool("save_reference_turntable") && 
          res_blk->get_block("save_reference_turntable_hydra_settings"))
        render_model_turntable(*(res_blk->get_block("save_reference_turntable_hydra_settings")), reference.model);

      if (res_blk->get_bool("save_model"))
      {
        assert(reconstructed_model.is_valid());
        assert(reconstructed_model.models.size() == 1); //TODO: save complex models with different materials too
        model_loader::save_model_to_obj(reconstructed_model.models[0], save_directory + "reconstructed_model.obj");

        if (reference.model.is_valid())
        {
          assert(reference.model.models.size() == 1); //TODO: save complex models with different materials too
          model_loader::save_model_to_obj(reference.model.models[0], save_directory + "reference_model.obj");          
        }
      } 
    }

    //print results of the reconstruction
    if (opt_blk->get_bool("verbose"))
    {
      debug("UPG Reconstruction finished\n");
      for (int i=0;i<opt_res.size();i++)
      {
        debug("=============================\n");
        debug("Optimization result %d\n",i);
        debug("Optimizer's loss    :%7.5f\n", opt_res[i].loss_optimizer);
        debug("Image-based quality :%7.5f\n", opt_res[i].quality_ir);
        if (reference.model.is_valid())
        debug("Model-based quality :%7.5f\n", opt_res[i].quality_synt);
        else
        debug("Model-based quality : ----- \n");
      }
      debug("=============================\n");
    }
    
    return opt_res;
  }
}