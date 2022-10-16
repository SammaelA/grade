#include "diff_optimization.h"
#include "diff_geometry_generation.h"
#include "mitsuba_python_interaction.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"
#include <functional>
#include <chrono>
#include <algorithm>
#include "common_utils/blk.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/silhouette.h"

namespace dopt
{
  class Optimizer
  {
  public:
    Optimizer(){};
    virtual ~Optimizer(){};
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) = 0;
  };

  class GradientDescentSimple : public Optimizer
  {
  public:
    GradientDescentSimple(float _alpha = 0.1) 
    {
      assert(_alpha > 0);
      alpha = _alpha;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
    {
      assert(x_prev.size() == x_grad.size());
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        x[i] = x_prev[i] - alpha*x_grad[i];
      }

      return x;
    }
  private:
    float alpha = 1;
  };

  class RMSprop : public Optimizer
  {
  public:
    RMSprop(float _alpha = 0.01, float _beta = 0.99, float _eps = 1e-8)
    {
      assert(_alpha > 0);
      assert(_beta > 0);
      assert(_beta < 1);
      assert(_eps > 0);

      alpha = _alpha;
      beta = _beta;
      eps = _eps;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
    {
      assert(x_prev.size() == x_grad.size());
      if (S.empty())
        S = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == S.size());
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        S[i] = beta * S[i] + (1-beta)*x_grad[i]*x_grad[i];
        x[i] = x_prev[i] - alpha*x_grad[i]/(sqrt(S[i]) + eps);
      }

      return x;
    }
  private:
    std::vector<float> S; 
    float alpha = 1;
    float beta = 1;
    float eps = 1;
  };

  class Adam : public Optimizer
  {
  public:
    Adam(float _alpha = 0.01, float _beta_1 = 0.9, float _beta_2 = 0.999, float _eps = 1e-8)
    {
      assert(_alpha > 0);
      assert(_beta_1 > 0);
      assert(_beta_1 < 1);
      assert(_beta_2 > 0);
      assert(_beta_2 < 1);
      assert(_eps > 0);
      
      alpha = _alpha;
      beta_1 = _beta_1;
      beta_2 = _beta_2;
      eps = _eps;
    }
    virtual std::vector<float> step(const std::vector<float> &x_prev, const std::vector<float> &x_grad, float value) override
    {
      assert(x_prev.size() == x_grad.size());
      if (S.empty())
        S = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == S.size());
      if (V.empty())
        V = std::vector<float>(x_prev.size(), 0);
      else
        assert(x_prev.size() == V.size());

      iter++;
      std::vector<float> x(x_prev.size(),0);
      for (int i=0;i<x_prev.size();i++)
      {
        V[i] = beta_1 * V[i] + (1-beta_1)*x_grad[i];
        float Vh = V[i] / (1 - pow(beta_1, iter)); 
        S[i] = beta_2 * S[i] + (1-beta_2)*x_grad[i]*x_grad[i];
        float Sh = S[i] / (1 - pow(beta_2, iter)); 

        x[i] = x_prev[i] - alpha*Vh/(sqrt(Sh) + eps);
      }

      return x;
    }
  private:
    std::vector<float> V; 
    std::vector<float> S; 
    float alpha = 1;
    float beta_1 = 1;
    float beta_2 = 1;
    float eps = 1;
    int iter = 0;
  };

  class Adam2 : public Optimizer
  {
  public:
    Adam2(float _alpha = 0.01, std::vector<float> _mask = {}, float _beta_1 = 0.9, float _beta_2 = 0.999, float _eps = 1e-8)
    {
      assert(_alpha > 0);
      assert(_beta_1 > 0);
      assert(_beta_1 < 1);
      assert(_beta_2 > 0);
      assert(_beta_2 < 1);
      assert(_eps > 0);
      
      mask = _mask;
      alpha = _alpha;
      beta_1 = _beta_1;
      beta_2 = _beta_2;
      eps = _eps;
    }
    virtual std::vector<float> step(const std::vector<float> &cur_x, const std::vector<float> &x_grad, float value) override
    {
      assert(cur_x.size() == x_grad.size());
      if (S.empty())
        S = std::vector<float>(cur_x.size(), 0);
      else
        assert(cur_x.size() == S.size());
      if (V.empty())
        V = std::vector<float>(cur_x.size(), 0);
      else
        assert(cur_x.size() == V.size());
      if (mask.empty())
        mask = std::vector<float>(cur_x.size(), 1);
      else
        assert(cur_x.size() == mask.size());

      std::vector<float> next_x(cur_x.size(),0);

      if (value < 1.2*prev_val)
      {
        iter++;
        for (int i=0;i<cur_x.size();i++)
        {
          V[i] = beta_1 * V[i] + (1-beta_1)*x_grad[i];
          float Vh = V[i] / (1 - pow(beta_1, iter)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*x_grad[i]*x_grad[i];
          float Sh = S[i] / (1 - pow(beta_2, iter)); 

          next_x[i] = cur_x[i] - mask[i]*alpha*Vh/(sqrt(Sh) + eps);
        }
        prev_x = cur_x; 
        prev_val = value;
      }
      else
      {
        for (int i=0;i<cur_x.size();i++)
        {
          float rnd = mask[i]*urand();
          next_x[i] = (1-rnd)*prev_x[i] + rnd*cur_x[i];
        } 
      }
      return next_x;
    }
  private:
    std::vector<float> V; 
    std::vector<float> S; 
    std::vector<float> prev_x;
    std::vector<float> mask;
    float prev_val = 1e9;
    float alpha = 1;
    float beta_1 = 1;
    float beta_2 = 1;
    float eps = 1;
    int iter = 0;
  };

  class DiffFunctionEvaluator
  {
  public:
    ~DiffFunctionEvaluator()
    {
      for (auto f : functions)
      {
        if (f)
          delete f;
      }
    }
    void init(dgen::generator_func _model_creator, int _model_creator_params_cnt)
    {
      model_creator = _model_creator;
      model_creator_params_cnt = _model_creator_params_cnt;
    }
    std::vector<float> get(const std::vector<float> &params)
    {
      return functions[find_or_add(params)]->Forward(0, params); 
    }
    std::vector<float> get_jac(const std::vector<float> &params)
    {
      return functions[find_or_add(params)]->Jacobian(params); 
    }
  private:
    int find_or_add(const std::vector<float> &params)
    {
      int hash = get_function_hash(params);
      auto it = hash_to_function_pos.find(hash);
      if (it == hash_to_function_pos.end())
      {
        std::vector<dgen::dfloat> X(params.size());
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        model_creator(X, Y);
        dgen::transform_by_scene_parameters(X, model_creator_params_cnt, Y);
        CppAD::ADFun<float> *f = new CppAD::ADFun<float>(X, Y); 
        functions.push_back(f);
        output_sizes.push_back(Y.size());
        int f_pos = functions.size()-1;
        hash_to_function_pos.emplace(hash, f_pos);

        return f_pos;
      }
      return it->second;
    }

    int get_function_hash(const std::vector<float> &params)
    {
      return 0;
    }
    int model_creator_params_cnt = 0;
    dgen::generator_func model_creator;
    std::map<int, int> hash_to_function_pos;
    std::vector<CppAD::ADFun<float> *> functions;
    std::vector<int> output_sizes;//same size as functions vector
  };

  struct OptimizationUnitGD
  {
    void init(int _id, const std::vector<float> &init_params, DiffFunctionEvaluator &_func, MitsubaInterface &_mi, 
              const std::vector<float> &_params_min, const std::vector<float> &_params_max,
              const std::vector<float> params_mask = {},
              bool _verbose = false)
    {
      assert(init_params.size() > 0);
      assert(init_params.size() == _params_min.size());
      assert(init_params.size() == _params_max.size());
      params_min = _params_min;
      params_max = _params_max;
      id = _id;
      verbose = _verbose;
      func = &_func;
      mi = &_mi;
      params = init_params;
      opt = new Adam2(0.015, params_mask);
      x_n = init_params.size();
      for (int i=0;i<8;i++)
        timers[i] = 0;
    }
    void iterate()
    {
      if (!func || !mi)
        return;
      if (verbose)
      {
        debug("params [");
        for (int j=0;j<x_n;j++)
        {
          debug("%.3f, ", params[j]);
        }
        debug("]\n");
      }
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    //float ms = 1e-4 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      std::vector<float> jac = func->get_jac(params);
      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      std::vector<float> res = func->get(params); 
      std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      float loss = mi->render_and_compare(res);
      std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
      mi->compute_final_grad(jac, x_n, res.size()/FLOAT_PER_VERTEX, final_grad);
      std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
      if (verbose)
      {
        debug("[%d] loss = %.3f, quality = %.3f\n", iterations, loss, get_quality());

        debug("grad {");
        for (int j=0;j<x_n;j++)
        {
          debug("%.3f ", final_grad[j]);
        }
        debug("}\n");
      }
      std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();

      if (loss < best_error)
      {
        best_error = loss;
        best_params = params;
        best_error_iter = iterations;
      }
      best_error_stat.push_back(best_error);
      quality = get_quality();
      params = opt->step(params, final_grad, loss);
      for (int i=0;i<params.size();i++)
      {
        if (params[i] <= params_min[i])
        {
          float rnd = urand(0, 1);
          params[i] = params_min[i] + rnd*MIN(params_min[i] - params[i], params_max[i] - params_min[i]);
        }
        else if (params[i] >= params_max[i])
        {
          float rnd = urand(0, 1);
          params[i] = params_max[i] - rnd*MIN(params[i] - params_max[i], params_max[i] - params_min[i]);
        }
      }
      iterations++;

      std::chrono::steady_clock::time_point t8 = std::chrono::steady_clock::now();
      timers[0] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      timers[1] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
      timers[2] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t5 - t3).count();
      timers[3] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
      timers[4] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
      timers[5] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
    }

    ~OptimizationUnitGD()
    {
      if (opt)
        delete opt;
    }

    float get_quality()
    {
      float time_stat = (iterations <= 10) ? 0.2*(15 - iterations) : (1.1 - 0.01*iterations);
      float val_stat = 1/(0.001 + best_error);
      float change_stat = 0;
      if (iterations >= 10)
        change_stat = 1 - best_error_stat[best_error_stat.size()-1 - 10] / (1e-9 + best_error);
      float decay_stat = MAX(0, 0.025*(iterations - best_error_iter - 10));
      return pow(val_stat, 2)*exp(time_stat)*exp(change_stat)*exp(-decay_stat);
    }

    void print_current_state()
    {
      debug("[%d][%d] q = %.3f, l = %.3f\n", id, iterations, quality, best_error);
    }

    void print_stat()
    {
      float total_time = 1e-6;
      for (int i=0;i<8;i++)
        total_time += timers[i];
      total_time *= 1e-3;

      debug("Optimization Unit %d statistics\n", id);
      debug("Input parameters count: %d\n", x_n);
      debug("Itarations: %d\n", iterations);
      debug("Best value: %.4f\n", best_error);
      debug("Total time: %.3f s\n", total_time);
      std::vector<std::string> markers = {
        "Model jacobian calculation", "Model calculation", "Rendering", "Final Gradient calculation",
        "Debug print", "Optimizaiton"
      };
      for (int i=0;i<markers.size();i++)
      {
        float time = 1e-3*timers[i];
        float part = 100*time/total_time;
        float per_iter = 1000*time/iterations;
        debug("%s: %.3f s (%.1f %), %.1f ms per iteration\n", markers[i].c_str(), time, part, per_iter);
      }
    }

    DiffFunctionEvaluator *func = nullptr;
    MitsubaInterface *mi = nullptr;
    Optimizer *opt = nullptr;
    std::vector<float> params;
    std::vector<float> best_params;
    std::vector<float> params_min;
    std::vector<float> params_max;
    float best_error = 1;
    std::vector<float> best_error_stat;
    int best_error_iter = 0;
    int iterations = 0;
    int x_n = 0;
    bool verbose = false;
    double timers[8];
    int id = 0;
    float quality = 10000;
  };

  class UShortVecComparator
  {
  public:
    bool operator()(const std::vector<unsigned short> &v1, const std::vector<unsigned short> &v2)
    {
      for (int i = 0; i < MIN(v1.size(), v2.size()); i++)
      {
        if (v1[i] < v2[i])
          return true;
        else if (v1[i] > v2[i])
          return false;
      }
      return false;
    }
  };

  void test()
  {
    Block gen_params, scene_params;
    std::vector<float> params_min, params_max;
    std::vector<unsigned short> init_bins_count;
    std::vector<unsigned short> init_bins_positions;
    load_block_from_file("dishes_gen_parameters_description.blk", gen_params);
    load_block_from_file("diff_gen_scene_parameters_description.blk", scene_params);

    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();
    int have_init_bins_cnt = 0;
    auto process_blk = [&](Block &blk){
      for (int i=0;i<blk.size();i++)
      {
        Block *pb = blk.get_block(i);
        if (!pb && pb->size()>0)
        {
          logerr("invalid parameter description\"%s\". It should be non-empty block", blk.get_name(i));
        }
        else
        {
          glm::vec2 min_max = pb->get_vec2("values", glm::vec2(1e9,-1e9));
          if (min_max.x > min_max.y)
            logerr("invalid parameter description\"%s\". It should have values:p2 with min and max values", blk.get_name(i));
          params_min.push_back(min_max.x);
          params_max.push_back(min_max.y);
          int bins_cnt = pb->get_int("init_bins_count", 0);
          if (bins_cnt<0 || bins_cnt>512)
          {
            bins_cnt = 0;
            logerr("invalid parameter description\"%s\". Bin count should be in [0, 512] interval", blk.get_name(i));
          }
          init_bins_count.push_back((unsigned short)bins_cnt);
          if (bins_cnt>0)
          {
            init_bins_positions.push_back(init_bins_count.size()-1);
            have_init_bins_cnt++;
          }
        }
      }
    };
    process_blk(gen_params);
    process_blk(scene_params);

    size_t x_n = gen_params_cnt + scene_params_cnt;
    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene). %d need start point selection\n", 
          x_n, gen_params_cnt, scene_params_cnt, have_init_bins_cnt);
    
    std::vector<float> reference_params{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                                        0.05, 0.35, 0.35, //hand params
                                        PI/5, PI, 0, 0, 0, 0};//rotation and transform
    //reference_params = std::vector<float>{4.281, 4.277, 4.641, 4.702, 4.639, 4.102, 3.748, 3.413, 3.771, 0.079, 0.013, 0.051, 0.659, 2.917, 0.098, 0.192, 0.210, 0.054};
    std::vector<float> init_params{4, 4, 4, 4, 4, 4, 4, 4, 4,
                                   0.05, 0.3, 0.4,
                                   0, PI, 0, 0, 0, 0};
    std::vector<float> params_mask{1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1,
                                   1, 1, 1, 1, 1, 1};
    int ref_image_size = 1024;
    int sel_image_size = 196;
    DiffFunctionEvaluator func;
    func.init(dgen::create_cup, gen_params_cnt);

    std::vector<float> reference = func.get(reference_params);

    MitsubaInterface mi("scripts", "emb_test");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(reference, "saves/reference.png");
    Texture t = engine::textureManager->load_unnamed_tex("saves/reference.png");
    SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
    Texture tex = se.get_silhouette(t, sel_image_size, sel_image_size);
    engine::textureManager->save_png_directly(tex, "saves/reference.png");

    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(sel_image_size, sel_image_size, 1, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE));
    mi.init_optimization("saves/reference.png", MitsubaInterface::LOSS_MSE, 1 << 16, false);
    /*
    Texture t = engine::textureManager->load_unnamed_tex("resources/textures/cup2.jpg");
    SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
    Texture tex = se.get_silhouette(t, mi.render_settings.image_w, mi.render_settings.image_h);
    engine::textureManager->save_png_directly(tex, "saves/reference.png");
    */
    bool simple_search = true;
    std::vector<float> best_params = init_params;
    float best_err = 1000;
    int total_iters = 0;

    if (simple_search)
    {
      OptimizationUnitGD opt_unit;
      opt_unit.init(0, init_params, func, mi, params_min, params_max, params_mask);
      for (int j=0;j<10;j++)
      {
          opt_unit.iterate();
          opt_unit.print_current_state();
      }
      opt_unit.print_stat();
      best_params = opt_unit.best_params;
      best_err = opt_unit.best_error;
      total_iters = opt_unit.iterations;
    }
    else
    {
    int full_cnt = 4;
    int use_cnt = 8;
    int base_iters = 500;
    int gd_iters = 5;
    std::map<std::vector<unsigned short>, int, UShortVecComparator> opt_unit_by_init_value_bins;
    std::vector<OptimizationUnitGD> opt_units(full_cnt);
    std::vector<std::vector<unsigned short>> opt_unit_bins(full_cnt);
    std::vector<int> indices_to_sort(full_cnt);
    int next_unit_id = 0;
    for (int i=0;i<full_cnt;i++)
    {
      indices_to_sort[i] = i;
      std::vector<unsigned short> descr(have_init_bins_cnt, 0);
      bool searching = true;
      int tries = 0;
      while (searching && tries<1000)
      {
        for (int j=0;j<have_init_bins_cnt;j++)
        {
          int max_bins = init_bins_count[init_bins_positions[j]];
          int bin = urandi(0, max_bins);
          descr[j] = bin;
        }
        tries++;
        if (opt_unit_by_init_value_bins.find(descr) == opt_unit_by_init_value_bins.end())
          searching = false;
      }

      if (!searching)
      {
        opt_unit_by_init_value_bins.emplace(descr, next_unit_id);
        std::vector<float> params = init_params;
        for (int j=0;j<have_init_bins_cnt;j++)
        {
          int pos = init_bins_positions[j];
          float val_from = params_min[pos] + descr[j]*(params_max[pos] - params_min[pos])/init_bins_count[pos];
          float val_to = params_min[pos] + (descr[j]+1)*(params_max[pos] - params_min[pos])/init_bins_count[pos];
          params[pos] = urand(val_from, val_to);
        }
        opt_unit_bins[i] = descr;
        opt_units[i].init(next_unit_id, params, func, mi, params_min, params_max, params_mask);
        next_unit_id++;
      }
    }
    for (int i=0;i<base_iters;i++)
    {
      /*choose a few best
      int iter_cnt = MAX(2, use_cnt - 0.33*i);
      if (i < 3)
        iter_cnt = full_cnt;
      for (int j=0;j<iter_cnt;j++)
      {
        for (int k=0;k<gd_iters;k++)
          opt_units[indices_to_sort[j]].iterate();
      }
      */
     //choose random, chance proportional to quality
      std::vector<double> sums;
      for (auto &unit : opt_units)
      {
        if (sums.empty())
          sums.push_back(unit.quality);
        else
          sums.push_back(sums.back() + unit.quality);
      }
      double rnd = urand(0, sums.back());
      for (int j=0;j<sums.size();j++)
      {
        if (sums[j] > rnd)
        {
          for (int k=0;k<gd_iters;k++)
            opt_units[j].iterate();
          break;
        }
      }

      std::sort(indices_to_sort.begin(), indices_to_sort.end(), 
                [&](const int& a, const int& b) -> bool{return opt_units[a].quality > opt_units[b].quality;});
      for (int ind : indices_to_sort)
      {
        opt_units[ind].print_current_state();
      }
      debugnl();
    }
    for (auto &unit : opt_units)
    {
      total_iters += unit.iterations;
      if (unit.best_error < best_err)
      {
        best_err = unit.best_error;
        best_params = unit.best_params;
      }
    }
    }
    
    std::vector<float> best_model = func.get(best_params);
    mi.render_model_to_file(best_model, "saves/selected.png");
    debug("Model optimization finished. %d iterations total. Best result saved to \"saves/selected.png\"\n", total_iters);
    debug("Best error: %f\n", best_err);
    debug("Best params: [");
    for (int j = 0; j < x_n; j++)
    {
      debug("%.3f, ", best_params[j]);
    }
    debug("]\n");
    mi.finish();
  }
}