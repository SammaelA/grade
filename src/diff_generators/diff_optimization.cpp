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
#include "save_utils/csv.h"
#include "graphics_utils/model_texture_creator.h"
#include "graphics_utils/modeling.h"

namespace dopt
{
  class UShortVecComparator
  {
  public:
    bool operator()(const std::vector<unsigned short> &v1, const std::vector<unsigned short> &v2) const
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
          next_x[i] = rnd*prev_x[i] + (1-rnd)*cur_x[i];
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
    void init(dgen::generator_func _model_creator, int _model_creator_params_cnt, std::vector<unsigned short> &variant_positions)
    {
      model_creator = _model_creator;
      model_creator_params_cnt = _model_creator_params_cnt;
      variant_params_positions = variant_positions;
    }
    std::vector<float> get(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      return functions[find_or_add(params, mq)]->Forward(0, params); 
    }
    std::vector<float> get_jac(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      return functions[find_or_add(params, mq)]->Jacobian(params); 
    }
  private:
    int find_or_add(const std::vector<float> &params, dgen::ModelQuality mq)
    {
      auto vs = get_variant_set(params);
      vs.push_back((unsigned short)mq.create_only_position);
      vs.push_back((unsigned short)mq.quality_level);
      auto it = variant_set_to_function_pos.find(vs);
      if (it == variant_set_to_function_pos.end())
      {
        debug("added new function {");
        for (auto &v : vs)
          debug("%d ", (int)v);
        debug("}\n");
        std::vector<dgen::dfloat> X(params.size());
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        model_creator(X, Y, mq);
        dgen::transform_by_scene_parameters(X, model_creator_params_cnt, Y);
        CppAD::ADFun<float> *f = new CppAD::ADFun<float>(X, Y); 
        functions.push_back(f);
        output_sizes.push_back(Y.size());
        int f_pos = functions.size()-1;
        variant_set_to_function_pos.emplace(vs, f_pos);

        return f_pos;
      }
      return it->second;
    }

    std::vector<unsigned short> get_variant_set(const std::vector<float> &params)
    {
      std::vector<unsigned short> vs;
      for (auto &pos : variant_params_positions)
      {
        vs.push_back((unsigned short)round(params[pos]));
      }
      return vs;
    }
    int model_creator_params_cnt = 0;
    dgen::generator_func model_creator;
    std::map<std::vector<unsigned short>, int, UShortVecComparator> variant_set_to_function_pos;
    std::vector<CppAD::ADFun<float> *> functions;
    std::vector<int> output_sizes;//same size as functions vector
    std::vector<unsigned short> variant_params_positions;
  };

  struct OptimizationUnitGD
  {
    void init(int _id, const std::vector<float> &init_params, DiffFunctionEvaluator &_func, MitsubaInterface &_mi, 
              const std::vector<float> &_params_min, const std::vector<float> &_params_max,
              CppAD::ADFun<float> *pr = nullptr,
              const std::vector<float> params_mask = {},
              bool _verbose = false)
    {
      assert(init_params.size() > 0);
      assert(init_params.size() == _params_min.size());
      assert(init_params.size() == _params_max.size());
      params_min = _params_min;
      params_max = _params_max;
      params_regularizer = pr;
      id = _id;
      verbose = _verbose;
      func = &_func;
      mi = &_mi;
      params = init_params;
      best_params = init_params;
      opt = new Adam2(0.05, params_mask);
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
      std::vector<float> jac = func->get_jac(params, dgen::ModelQuality(true, 0));
      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      std::vector<float> res = func->get(params, dgen::ModelQuality(true, 0)); 
      std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
      std::vector<float> final_grad = std::vector<float>(x_n, 0);
      float loss = mi->render_and_compare(res, timers);
      std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
      mi->compute_final_grad(jac, x_n, res.size()/FLOAT_PER_VERTEX, final_grad);
      std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
      float reg_q = 0.33;
      if (params_regularizer)
      {
        std::vector<float> reg_res = params_regularizer->Forward(0, params);
        std::vector<float> reg_jac = params_regularizer->Jacobian(params);
        loss += reg_q*reg_res[0];
        for (int i=0;i<MIN(final_grad.size(), reg_jac.size());i++)
          final_grad[i] += reg_q*reg_jac[i];
      }
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
      //timers[2] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t5 - t3).count();
      timers[5] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
      timers[6] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
      timers[7] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
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
      debug("Total time: %.3f s, %.1f ms per iteration\n", total_time, 1000*total_time/iterations);
      std::vector<std::string> markers = {
        "Model jacobian calculation", "Model calculation", "Send model to mitsuba", "Rendering", "Get result from mitsuba",
        "Final Gradient calculation", "Debug print", "Optimizaiton"
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
    CppAD::ADFun<float> *params_regularizer = nullptr;
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

  struct OptimizationResult
  {
    std::vector<float> best_params;
    float best_err;
    int total_iters;
  };

  void optimizer_simple_search(Block *settings, CppAD::ADFun<float> *f_reg, DiffFunctionEvaluator &func, MitsubaInterface &mi,
                       std::vector<float> &init_params, std::vector<float> &params_min, std::vector<float> &params_max,
                       std::vector<float> &params_mask, int verbose_level, std::string save_stat_path,
                       std::vector<unsigned short> init_bins_count, std::vector<unsigned short> init_bins_positions,
                       std::vector<std::vector<float>> parameter_presets,
                       OptimizationResult &opt_result)
  {
    int iterations = settings->get_int("iterations", 40);
    OptimizationUnitGD opt_unit;
    opt_unit.init(0, init_params, func, mi, params_min, params_max, f_reg, params_mask, verbose_level == 2);
    for (int j = 0; j < iterations; j++)
    {
      opt_unit.iterate();
      if (verbose_level > 0)
        opt_unit.print_current_state();
    }
    if (verbose_level > 0)
      opt_unit.print_stat();

    opt_result.best_params = opt_unit.best_params;
    opt_result.best_err = opt_unit.best_error;
    opt_result.total_iters = opt_unit.iterations;

    if (save_stat_path != "")
    {
      // save optimization staticstics to csv
      CSVData csv = CSVData({"iteration", "unit_id", "loss"});
      int i = 0;
      for (auto &val : opt_unit.best_error_stat)
      {
        std::vector<float> row{(float)i, (float)opt_unit.id, val};
        csv.add_row(row);
        i++;
      }
      CSVSaver saver;
      saver.save_csv_in_file(csv, save_stat_path);
    }
  }

  std::vector<float> get_new_init_point(std::vector<float> &init_params, std::vector<float> &params_min, std::vector<float> &params_max,
                                        std::vector<unsigned short> init_bins_count, std::vector<unsigned short> init_bins_positions,
                                        std::map<std::vector<unsigned short>, int, UShortVecComparator> opt_unit_by_init_value_bins,
                                        int unit_id, bool fill_rest_with_random)
  {
    std::vector<unsigned short> descr(init_bins_count.size(), 0);
    bool searching = true;
    int tries = 0;
    while (searching && tries < 1000)
    {
      for (int j = 0; j < init_bins_count.size(); j++)
      {
        int max_bins = init_bins_count[j];
        int bin = urandi(0, max_bins);
        descr[j] = bin;
      }
      tries++;
      if (opt_unit_by_init_value_bins.find(descr) == opt_unit_by_init_value_bins.end())
        searching = false;
    }
    debug("%d descr [", (int)searching);
    for (auto &d : descr)
      debug("%d ", (int)d);
    debug("]\n");
    opt_unit_by_init_value_bins.emplace(descr, unit_id);
    std::vector<float> params = init_params;
    if (fill_rest_with_random)
    {
      for (int i=0; i<params.size();i++)
      {
        params[i] = params_min[i] + ((urand(0.25, 0.75)+urand(0.25, 0.75))/2)*(params_max[i] - params_min[i]);
      }
    }
    for (int j = 0; j < init_bins_count.size(); j++)
    {
      int pos = init_bins_positions[j];
      float val_from = params_min[pos] + descr[j] * (params_max[pos] - params_min[pos]) / init_bins_count[j];
      float val_to = params_min[pos] + (descr[j] + 1) * (params_max[pos] - params_min[pos]) / init_bins_count[j];
      params[pos] = urand(val_from, val_to);
    }

    return params;
  }

  void optimizer_advanced_search(Block *settings, CppAD::ADFun<float> *f_reg, DiffFunctionEvaluator &func, MitsubaInterface &mi,
                        std::vector<float> &init_params, std::vector<float> &params_min, std::vector<float> &params_max,
                        std::vector<float> &params_mask, int verbose_level, std::string save_stat_path,
                        std::vector<unsigned short> init_bins_count, std::vector<unsigned short> init_bins_positions,
                        std::vector<std::vector<float>> parameter_presets,
                        OptimizationResult &opt_result)
  {
    int full_cnt = settings->get_int("start_points", 4);
    int base_iters = settings->get_int("base_iterations", 200);
    int gd_iters = settings->get_int("gd_iterations_per_base_iteration", 1);

    std::map<std::vector<unsigned short>, int, UShortVecComparator> opt_unit_by_init_value_bins;
    std::vector<OptimizationUnitGD> opt_units(full_cnt);
    std::vector<int> indices_to_sort(full_cnt);
    int next_unit_id = 0;
    for (int i = 0; i < full_cnt; i++)
    {
      std::vector<float> cur_init_params = init_params;
      if (parameter_presets.size() > 0 && urand() > 0.5)
      {
        cur_init_params = parameter_presets[(int)urandi(0, parameter_presets.size())];
      }
      std::vector<float> params = get_new_init_point(cur_init_params, params_min, params_max, init_bins_count, init_bins_positions, 
                                                     opt_unit_by_init_value_bins, next_unit_id, false);
      indices_to_sort[i] = i;
      opt_units[i].init(next_unit_id, params, func, mi, params_min, params_max, f_reg, params_mask, verbose_level == 2);
      next_unit_id++;
    }
    for (int i = 0; i < base_iters; i++)
    {
      /*choose a few best
      int use_cnt = 8;
      int iter_cnt = MAX(2, use_cnt - 0.33*i);
      if (i < 3)
        iter_cnt = full_cnt;
      for (int j=0;j<iter_cnt;j++)
      {
        for (int k=0;k<gd_iters;k++)
          opt_units[indices_to_sort[j]].iterate();
      }
      */
      // choose random, chance proportional to quality
      std::vector<double> sums;
      for (auto &unit : opt_units)
      {
        if (sums.empty())
          sums.push_back(unit.quality);
        else
          sums.push_back(sums.back() + unit.quality);
      }
      double rnd = urand(0, sums.back());
      for (int j = 0; j < sums.size(); j++)
      {
        if (sums[j] > rnd)
        {
          for (int k = 0; k < gd_iters; k++)
            opt_units[j].iterate();
          break;
        }
      }

      std::sort(indices_to_sort.begin(), indices_to_sort.end(),
                [&](const int &a, const int &b) -> bool
                { return opt_units[a].quality > opt_units[b].quality; });
      if (verbose_level > 0)
      {
        for (int ind : indices_to_sort)
        {
          opt_units[ind].print_current_state();
        }
        debugnl();
      }
    }
    for (auto &unit : opt_units)
    {
      opt_result.total_iters += unit.iterations;
      if (unit.best_error < opt_result.best_err)
      {
        opt_result.best_err = unit.best_error;
        opt_result.best_params = unit.best_params;
      }
      if (verbose_level > 0)
        unit.print_stat();
    }

    if (save_stat_path != "")
    {
      // save optimization staticstics to csv
      CSVData csv = CSVData({"iteration", "unit_id", "loss"});
      for (auto &unit : opt_units)
      {
        int i = 0;
        for (auto &val : unit.best_error_stat)
        {
          std::vector<float> row{(float)i, (float)unit.id, val};
          csv.add_row(row);
          i++;
        }
      }
      CSVSaver saver;
      saver.save_csv_in_file(csv, save_stat_path);
    }
  }

  float get_quality_for_memetic(OptimizationUnitGD &unit, int ga_iter, float decay_q = 0.005, float iter_q = 0.2)
  {
    float val_stat = 1 / (0.001 + unit.best_error);
    float decay_stat = MAX(0, decay_q*(unit.iterations - unit.best_error_iter));
    return pow(val_stat, 1+iter_q*ga_iter) * exp(-decay_stat);
  }

  void optimizer_memetic(Block *settings, CppAD::ADFun<float> *f_reg, DiffFunctionEvaluator &func, MitsubaInterface &mi,
                         std::vector<float> &init_params, std::vector<float> &params_min, std::vector<float> &params_max,
                         std::vector<float> &params_mask, int verbose_level, std::string save_stat_path,
                         std::vector<unsigned short> init_bins_count, std::vector<unsigned short> init_bins_positions,
                         std::vector<std::vector<float>> parameter_presets,
                         OptimizationResult &opt_result)
  {
    int ga_iters = settings->get_int("genetic_algorithm_iterations", 10);
    int gd_iters = settings->get_int("gradient_descent_iterations", 40);
    int initial_population_size = settings->get_int("initial_population_size", 8);
    float remove_chance = settings->get_double("remove_chance", 0.5);
    float recombination_chance = settings->get_double("recombination_chance", 0.2);
    float mutation_chance = settings->get_double("mutation_chance", 0.2);
    float reinit_chance = settings->get_double("reinit_chance", 0.0);
    float mutation_power = settings->get_double("mutation_power", 0.2);

    int next_unit_id = 0;
    int x_n = params_min.size();
    std::vector<OptimizationUnitGD> population(initial_population_size);
    std::map<std::vector<unsigned short>, int, UShortVecComparator> opt_unit_by_init_value_bins;

    auto init_random = [&](OptimizationUnitGD &unit)
    {
      std::vector<float> params = get_new_init_point(init_params, params_min, params_max, init_bins_count, init_bins_positions, 
                                                     opt_unit_by_init_value_bins, next_unit_id, false);
      unit.init(next_unit_id, params, func, mi, params_min, params_max, f_reg, params_mask, verbose_level == 2);
      next_unit_id++;
    };

    auto mutate = [&](OptimizationUnitGD &unit, OptimizationUnitGD &u1)
    {
      std::vector<float> params(params_max.size());
      params = u1.best_params;
      int pos = (int)(urandi(0, params_max.size()));
      params[pos] = params_min[pos] + ((urand(0.25, 0.75)+urand(0.25, 0.75))/2)*(params_max[pos] - params_min[pos]);
      /*
      for (int i=0; i<x_n;i++)
      {
        if (urand() < mutation_power)
          params[i] = urand(params_min[i], params_max[i]);
        else
          params[i] = u1.best_params[i];
      }
      */
      unit.init(next_unit_id, params, func, mi, params_min, params_max, f_reg, params_mask, verbose_level == 2);
      next_unit_id++;
    };

    auto recombime = [&](OptimizationUnitGD &unit, OptimizationUnitGD &u1, OptimizationUnitGD &u2)
    {
      std::vector<float> params(params_max.size());
      for (int i=0; i<x_n;i++)
      {
        if (urand() < 0.5)
          params[i] = u1.best_params[i];
        else
          params[i] = u2.best_params[i];
      }
      unit.init(next_unit_id, params, func, mi, params_min, params_max, f_reg, params_mask, verbose_level == 2);
      next_unit_id++;
    };

    //initialize population
    for (int i=0;i<initial_population_size;i++)
    {
      init_random(population[i]);
    }

    //main loop of GA
    for (int iter = 0; iter < ga_iters; iter++)
    {
      //improve current population with Gradient Descent
      float q_sum = 0;
      for (int i = 0; i < population.size(); i++)
        q_sum += population[i].id == -1 ? 0 : get_quality_for_memetic(population[i], iter);

      for (int unit_pos = 0; unit_pos < population.size(); unit_pos++)
      {
        auto &unit = population[unit_pos];
        if (unit.id != -1)
        {
          int iters = gd_iters * get_quality_for_memetic(unit, iter)/q_sum;
          if (unit.iterations == 0)
            iters = (float)gd_iters / population.size();
          for (int i=0;i<iters;i++)
            unit.iterate();
        }
      }

      //sort current population by quality
      std::vector<std::pair<int, float>> unit_indices(population.size());
      for (int i=0;i<population.size();i++)
      {
        float quality = population[i].id >= 0 ? get_quality_for_memetic(population[i], iter) : -1;
        unit_indices[i] = std::pair<int, float>(i, quality);
      }
      std::sort(unit_indices.begin(), unit_indices.end(), 
                [&](const std::pair<int, float>& a, const std::pair<int, float>& b) -> bool{return a.second > b.second;});
      
      if (verbose_level > 0)
      {
        debug("iteration %d\n", iter);
        for (int i=0;i<population.size();i++)
        {
          auto &unit = population[unit_indices[i].first];
          debug("[%d][%f][%f] - %f  ", unit.id, unit_indices[i].second, unit.best_error, unit.params[10]);
          for (int j=0;j<init_bins_count.size();j++)
          {
            int pos = init_bins_positions[j];
            int bin = init_bins_count[j]*(unit.params[pos] - params_min[pos])/(params_max[pos] - params_min[pos]);
            debug("%d ", bin);
          }
          debugnl();
        }
      }
      //if there are two or more units inside the same bin set, no need to save them both.
      //they will probably find the same local mimimum (with gradient descent). That's 
      //why bins are actually used. So we save only the best one

      std::vector<std::vector<unsigned short>> unit_bins;
      for (auto &unit : population)
      {
        std::vector<unsigned short> bins;
        for (int j=0;j<init_bins_count.size();j++)
        {
          int pos = init_bins_positions[j];
          int bin = init_bins_count[j]*(unit.params[pos] - params_min[pos])/(params_max[pos] - params_min[pos]);
          bins.push_back(bin);
        }
        unit_bins.push_back(bins);
      }
      for (int i=0;i<population.size();i++)
      {
        for (int j=i+1;j<population.size();j++)
        {
          if (population[i].id >= 0 && population[j].id >= 0)
          {
            bool same_bins = true;
            for (int b=0;b<init_bins_count.size();b++)
            {
              if (unit_bins[i][b] != unit_bins[j][b])
              {
                same_bins = false;
                break;
              }
            }
            if (same_bins && verbose_level)
            {
              debug("Units %d and %d have parameters in the same bins. Removing one of them\n", population[i].id, population[j].id);
              if (population[i].best_error < population[j].best_error)
                population[j].id = -1;
              else
                population[i].id = -1;
            }
          }
        }
      }
      //remove worst units and replace them with mutated, recombined or new random units
      int remove_cnt = MIN(population.size() - 2, round(remove_chance * population.size()));
      for (int i = 0; i < remove_cnt; i++)
      {
        //clear unit and mark it as invalid
        auto &u = population[population.size() - 1 - i];
        opt_result.total_iters += u.iterations;
        u = OptimizationUnitGD();
        u.id = -1;

        float rnd = urand();
        if (rnd < recombination_chance)
        {
          //choose two random units from existing ones and recombine them
          auto &u1 = population[(int)urandi(0, population.size() - remove_cnt)];
          auto &u2 = population[(int)urandi(0, population.size() - remove_cnt)];
          recombime(u, u1, u2);
        }
        else if (rnd < recombination_chance + mutation_chance)
        {
          auto &u1 = population[(int)urandi(0, population.size() - remove_cnt)];
          mutate(u, u1);
        }
        else if (rnd < recombination_chance + mutation_chance + reinit_chance)
        {
          init_random(u);
        }
      }
      for (auto &unit : population)
      {
        if (unit.id >= 0)
        {
          if (iter == ga_iters - 1)
            opt_result.total_iters += unit.iterations;
          if (unit.best_error < opt_result.best_err)
          {
            opt_result.best_err = unit.best_error;
            opt_result.best_params = unit.best_params;
          }
        }
      }
    }
  }

  void test()
  {
    std::vector<float> reference_params{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                                        0.8,// y_scale
                                        1, //has handle variant
                                        0.05, 0.35, 0.35, //hand params
                                        0, PI, 0, 0, 0, 0};//rotation and transform
    std::vector<float> init_params{4, 4, 4, 4, 4, 4, 4, 4, 4,
                                   1,
                                   1,
                                   0.05, 0.1, 0.1,
                                   0, PI, 0, 0, 0, 0};
    std::vector<float> params_mask{1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1,
                                   1,
                                   1, 1, 1,
                                   1, 1, 1, 1, 1, 1};
    Block settings_blk;
    settings_blk.add_string("parameters_description", "dishes_gen_parameters_description.blk");
    settings_blk.add_string("scene_description", "diff_gen_scene_parameters_description.blk");
    settings_blk.add_bool("synthetic_reference", true);
    settings_blk.add_string("reference_path", "resources/textures/cup1.jpg");
    settings_blk.add_arr("reference_params", reference_params);
    settings_blk.add_arr("init_params", init_params);
    settings_blk.add_arr("params_mask", params_mask);

    Block opt_settings_block;
    opt_settings_block.add_int("iterations", 100);
    settings_blk.add_block("simple_search", &opt_settings_block);

    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    image_based_optimization(settings_blk, mi);
  }

  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi)
  {
    Block gen_params, scene_params, presets_blk;
    std::vector<float> params_min, params_max;
    std::vector<unsigned short> init_bins_count;
    std::vector<unsigned short> init_bins_positions;
    std::vector<unsigned short> variant_count;
    std::vector<unsigned short> variant_positions;
    std::vector<std::vector<float>> parameter_presets;

    load_block_from_file(settings_blk.get_string("parameters_description"), gen_params);
    load_block_from_file(settings_blk.get_string("scene_description"), scene_params);
    load_block_from_file(settings_blk.get_string("presets_block"), presets_blk);
    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();
    size_t x_n = gen_params_cnt + scene_params_cnt;

    int verbose_level = settings_blk.get_int("verbose_level", 1);
    int ref_image_size = settings_blk.get_int("reference_image_size", 512);
    int sel_image_size = settings_blk.get_int("selection_image_size", 196);
    bool by_reference = settings_blk.get_bool("synthetic_reference", true);
    bool texture_extraction = settings_blk.get_bool("texture_extraction", false);
    std::string search_algorithm = settings_blk.get_string("search_algorithm", "simple_search");
    std::string save_stat_path = settings_blk.get_string("save_stat_path", "");
    std::string saved_result_path = settings_blk.get_string("saved_result_path", "saves/selected_final.png");
    std::string saved_textured_path = settings_blk.get_string("saved_textured_path", "saves/selected_textured.png");
    std::string saved_initial_path = settings_blk.get_string("saved_initial_path", "");
    std::string reference_path = settings_blk.get_string("reference_path", "");

    std::vector<float> reference_params, init_params, params_mask;
    if (by_reference)
    {
      settings_blk.get_arr("reference_params", reference_params);
      if (reference_params.size() != x_n)
      {
        logerr("DOpt Error: reference_params has %d values, it should have %d", reference_params.size(), x_n);
        return 1.0;
      }
    }

    settings_blk.get_arr("params_mask", params_mask);
    if (params_mask.empty())
      params_mask = std::vector<float>(x_n, 1);
    else if (params_mask.size() != x_n)
    {
      logerr("DOpt Error: params_mask has %d values, it should have %d", params_mask.size(), x_n);
      return 1.0;
    }

    settings_blk.get_arr("init_params", init_params);
    if (init_params.size() != x_n && init_params.size() > 0)
    {
      logerr("DOpt Error: init_params has %d values, it should have %d", params_mask.size(), x_n);
      return 1.0;
    }

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
          bool is_variant = pb->get_bool("is_variant", false);
          if (is_variant)
          {
            if (bins_cnt > 0)
            {
              logerr("invalid parameter description\"%s\". Variant parameter should not have explicit init_bins_count", blk.get_name(i));
            }
            int imin = round(min_max.x);
            int imax = round(min_max.y);
            if (abs((float)imin - min_max.x) > 1e-3 || abs((float)imax - min_max.y) > 1e-3)
            {
              logerr("invalid parameter description\"%s\". Variant parameter should have integer min max values", blk.get_name(i));
            }
            bins_cnt = imax - imin + 1;
            variant_count.push_back(bins_cnt);
            variant_positions.push_back(params_min.size()-1);
          }
          if (bins_cnt<0 || bins_cnt>512)
          {
            bins_cnt = 0;
            logerr("invalid parameter description\"%s\". Bin count should be in [0, 512] interval", blk.get_name(i));
          }

          if (bins_cnt>0)
          {
            init_bins_count.push_back((unsigned short)bins_cnt);
            init_bins_positions.push_back(params_min.size()-1);
          }
        }
      }
    };
    process_blk(gen_params);
    process_blk(scene_params);

    for (unsigned short &pos : variant_positions)
      params_mask[pos] = 0;

    if (init_params.empty())
    {
      init_params = std::vector<float>(x_n, 0);
      for (int i=0;i<x_n;i++)
      {
        init_params[i] = 0.5*(params_min[i] + params_max[i]);
      }
    }
    else 
    {
      for (int i=0;i<x_n;i++)
      {
        if (init_params[i] < params_min[i] || init_params[i] > params_max[i])
        {
          logerr("Wrong initial value %f for parameter %d. In should be in [%f, %f] interval",
                 init_params[i], i, params_min[i], params_max[i]);
          init_params[i] = CLAMP(init_params[i], params_min[i], params_max[i]);
        }
      }
    }

    if (!reference_params.empty())
    {
      for (int i=0;i<x_n;i++)
      {
        if (reference_params[i] < params_min[i] || reference_params[i] > params_max[i])
        {
          logerr("Wrong reference value %f for parameter %d. In should be in [%f, %f] interval",
                 reference_params[i], i, params_min[i], params_max[i]);
          reference_params[i] = CLAMP(reference_params[i], params_min[i], params_max[i]);
        }
      }
    }

    {
      //presets are manually created valid sets of generator parmeters, representing different types of objects,
      //that generator can create. They do not include scene parameters
      std::map<std::string, int> param_n_by_name;
      for (int i=0;i<gen_params.size();i++)
      {
        Block *pb = gen_params.get_block(i);
        if (pb)
        {
          param_n_by_name.emplace(gen_params.get_name(i), i);
        }
      }
      for (int i=0;i<presets_blk.size();i++)
      {
        Block *preset_block = presets_blk.get_block(i);
        if (preset_block)
        {
          std::vector<float> preset = init_params;
          if (preset_block->has_tag("compact"))
          {
            std::vector<float> params;
            preset_block->get_arr("params", params);
            if (params.size() == gen_params_cnt)
            {
              for (int k=0;k<gen_params_cnt;k++)
                preset[k] = params[k];
            }
            else
            {
              logerr("Error: compact preset %s has %d parameters, while the generator requests %d",
                     presets_blk.get_name(i).c_str(), params.size(), gen_params_cnt);
            }
          }
          else
          {
            for (int j=0;j<preset_block->size();j++)
            {
              auto it = param_n_by_name.find(preset_block->get_name(j));
              if (it != param_n_by_name.end())
              {
                float val = 0;
                if (preset_block->get_type(i) == Block::ValueType::DOUBLE)
                  val = preset_block->get_double(i);
                else if (preset_block->get_type(i) == Block::ValueType::INT)
                  val = preset_block->get_int(i);
                else if (preset_block->get_type(i) == Block::ValueType::BOOL)
                  val = (int)(preset_block->get_bool(i));
                else
                  logerr("parameter %s of preset %s has unknown parameter type. It should be double, int or bool",
                         preset_block->get_name(j).c_str(), presets_blk.get_name(i).c_str());
                preset[it->second] = val;
              }
              else
              {
                logerr("Unknown parameter name %s in preset %s",
                         preset_block->get_name(j).c_str(), presets_blk.get_name(i).c_str());
              }
            }
          }
          parameter_presets.push_back(preset);
          debug("read preset %s \n", presets_blk.get_name(i).c_str());
          for (int k=0;k<gen_params_cnt;k++)
            debug("%f ", preset[k]);
          debugnl();
        }
      }
    }

    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene). %d SP %d var\n", 
          x_n, gen_params_cnt, scene_params_cnt, init_bins_count.size(), variant_count.size());

    CppAD::ADFun<float> f_reg;
    {
      std::vector<dgen::dfloat> X(init_params.size());
      for (int i = 0; i < init_params.size(); i++)
        X[i] = init_params[i];
      std::vector<dgen::dfloat> Y;
      CppAD::Independent(X);
      Y.resize(1);
      Y[0] = dgen::parameters_limits_reg(X, params_min, params_max) + dgen::parameters_cup_reg(X);
      f_reg = CppAD::ADFun<float>(X, Y);
    }

    DiffFunctionEvaluator func;
    func.init(dgen::create_cup, gen_params_cnt, variant_positions);

    Texture reference_tex, reference_mask;

    if (by_reference)
    {
      std::vector<float> reference = func.get(reference_params);
      mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
      mi.render_model_to_file(reference, "saves/reference.png", dgen::ModelLayout());
      reference_tex = engine::textureManager->load_unnamed_tex("saves/reference.png");
      SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
      reference_mask = se.get_silhouette(reference_tex, sel_image_size, sel_image_size);
      engine::textureManager->save_png_directly(reference_mask, "saves/reference.png");
    }
    else
    {
      reference_tex = engine::textureManager->load_unnamed_tex(reference_path);
      SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
      reference_mask = se.get_silhouette(reference_tex, sel_image_size, sel_image_size);
      engine::textureManager->save_png_directly(reference_mask, "saves/reference.png");
    }
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(sel_image_size, sel_image_size, 1, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE));
    mi.init_optimization("saves/reference.png", MitsubaInterface::LOSS_MIXED, 1 << 16, dgen::ModelLayout(0, 3, 3, 3, 8), false);

    OptimizationResult opt_result{init_params, 1000, 0};

    Block *opt_settings = settings_blk.get_block(search_algorithm);
    if (!opt_settings)
    {
      logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
      return 1.0;
    }
    if (search_algorithm == "simple_search")
    {
      optimizer_simple_search(opt_settings, &f_reg, func, mi, init_params, params_min, params_max, params_mask, verbose_level, save_stat_path,
                              init_bins_count, init_bins_positions, parameter_presets, opt_result);
    }
    else if (search_algorithm == "advanced_search")
    {
      optimizer_advanced_search(opt_settings, &f_reg, func, mi, init_params, params_min, params_max, params_mask, verbose_level, save_stat_path,
                                init_bins_count, init_bins_positions, parameter_presets, opt_result);
    }
    else if (search_algorithm == "memetic")
    {
      optimizer_memetic(opt_settings, &f_reg, func, mi, init_params, params_min, params_max, params_mask, verbose_level, save_stat_path,
                        init_bins_count, init_bins_positions, parameter_presets, opt_result);
    }
    else
    {
      logerr("Unknown optimizer algorithm %s", search_algorithm.c_str());
      return 1.0;
    }
    
    std::vector<float> best_model = func.get(opt_result.best_params, dgen::ModelQuality(false, 3));
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(best_model, saved_result_path, dgen::ModelLayout());
    if (saved_initial_path != "")
    {
      std::vector<float> initial_model = func.get(init_params);
      mi.render_model_to_file(initial_model, saved_initial_path, dgen::ModelLayout());
    }
    if (texture_extraction)
    {
      ModelTex mt;
      Model *m = new Model();
      visualizer::simple_mesh_to_model_332(best_model, m);
      m->update();
      Texture res_tex = mt.getTexbyUV(reference_mask, *m, reference_tex, 3);
      engine::textureManager->save_png(res_tex, "reconstructed_tex");

      mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/reconstructed_tex.png"));
      mi.render_model_to_file(best_model, saved_textured_path, dgen::ModelLayout());
    }
    debug("Model optimization finished. %d iterations total. Best result saved to \"%s\"\n", opt_result.total_iters, saved_result_path.c_str());
    debug("Best error: %f\n", opt_result.best_err);
    debug("Best params: [");
    for (int j = 0; j < x_n; j++)
    {
      debug("%.3f, ", opt_result.best_params[j]);
    }
    debug("]\n");
    mi.finish();

    return opt_result.best_err;
  }
}