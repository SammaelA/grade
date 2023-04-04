#include "optimization.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"

namespace opt
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


  std::vector<float> get_new_init_point(std::vector<float> &init_params, std::vector<float> &params_min, std::vector<float> &params_max,
                                        std::vector<unsigned short> init_bins_count, std::vector<unsigned short> init_bins_positions,
                                        std::map<std::vector<unsigned short>, int, UShortVecComparator> &opt_unit_by_init_value_bins,
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

  void GridSearchAdam::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                                init_params_func &get_init_params)
  {
    std::vector<unsigned short> init_bins_count;
    std::vector<unsigned short> init_bins_positions;
    std::vector<std::vector<float>> init_bins_values;
    std::vector<std::string> params_names;

    settings.get_arr("params_names", params_names);
    Block *grid = settings.get_block("grid");
    if (!grid)
      logerr("GridSearchAdam: \"grid\" block not found");
    for (int i=0;i<grid->size();i++)
    {
      Block *axis = grid->get_block(i);
      if (axis)
      {
        int param_id = -1;
        for (int j=0;j<params_names.size();j++)
        {
          if (params_names[j] == grid->get_name(i))
          {
            param_id = j;
            break;
          }
        }
        if (param_id == -1)
          logerr("GridSearchAdam: grid has unknown parameter \"%s\"", grid->get_name(i).c_str());
        else
        {
          std::vector<float> values;
          axis->get_arr("values", values);
          if (values.size() > 0)
          {
            init_bins_positions.push_back(param_id);
            init_bins_count.push_back(values.size());
            init_bins_values.push_back(values);
          }
        }
      }
    }

    assert(min_X.size() > 0);
    assert(min_X.size() == max_X.size());
    assert(init_bins_count.size() == init_bins_positions.size());

    bool verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;
    int local_search_iterations = settings.get_int("local_search_iterations", 100);
    float local_search_learning_rate = settings.get_double("local_search_learning_rate", 0.25);
    int start_points_count = settings.get_int("start_points_count", 50);
    float grid_params_gradient_mult = settings.get_double("grid_params_gradient_mult", 1);

    std::vector<float> derivatives_mult(min_X.size(), 1);
    for (auto &i : init_bins_positions)
      derivatives_mult[i] = grid_params_gradient_mult;

    int N = min_X.size();

    if (init_bins_count.size() == 0)
    {
      debug("GridSearchAdam: init_bins_count vector is empty. Falling back to simple Adam\n");
      Optimizer *opt = new Adam();
      Block adam_settings;
      adam_settings.add_double("learning_rate", local_search_learning_rate);
      adam_settings.add_int("iterations", local_search_iterations*start_points_count);
      opt->optimize(F, min_X, max_X, adam_settings);
      best_params = opt->get_best_result(&best_result);

      delete opt;

      return;
    }
    std::vector<std::vector<unsigned short>> start_points_bins;
    uint64_t combinations = 1;
    for (auto &p : init_bins_count)
      combinations *= p;
    if (combinations <= start_points_count)
    {
      debug("GridSearchAdam: start_points_count (%d) is enough to test all combinations (%lu)\n",start_points_count, combinations);

      std::vector<unsigned short> combination = std::vector<unsigned short>(init_bins_count.size(), 0);
      for (int i=0;i<combinations;i++)
      {
        start_points_bins.push_back(combination);
        combination[0]++;
        
        int comb_pos = 0;
        while (combination[comb_pos] == init_bins_count[comb_pos])
        {
          for (int c = 0; c<=comb_pos; c++)
            combination[c] = 0;
          comb_pos++;
          if (comb_pos != init_bins_count.size())
            combination[comb_pos]++;
        }
      }
    }
    else
    {
      debug("GridSearchAdam: choosing %d starting points from %lu possible combinations\n",start_points_count, combinations);
      std::map<std::vector<unsigned short>, int, UShortVecComparator> opt_unit_by_init_value_bins;

      for (int i=0;i<start_points_count;i++)
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
        if (!searching)
          opt_unit_by_init_value_bins.emplace(descr, i);
        else
          logerr("GridSearchAdam: failed to find valid starting point in 1000 tries.");
      }

      for (auto &p : opt_unit_by_init_value_bins)
        start_points_bins.push_back(p.first);
    }

    int b = 0;
    for (auto &bins : start_points_bins)
    {
      debug("bin %d[", bins.size());
      for (auto &d : bins)
        debug("%d ", (int)d);
      debug("]\n");

      std::vector<float> start_params;
      settings.get_arr("initial_params", start_params);
      if (start_params.empty())
      {
        logerr("GSA: No initial params. Random values will be picked");
        start_params = std::vector<float>(N, 0);
        for (int i=0; i<start_params.size();i++)
        {
          start_params[i] = min_X[i] + ((urand(0, 1)+urand(0,1))/2)*(max_X[i] - min_X[i]);
        }
      }
      for (int j = 0; j < init_bins_count.size(); j++)
      {
        int pos = init_bins_positions[j];
        if (init_bins_values[j].size() > 0)
          start_params[pos] = init_bins_values[j][bins[j]];
        else
        {
          float val_from = min_X[pos] + bins[j] * (max_X[pos] - min_X[pos]) / init_bins_count[j];
          float val_to = min_X[pos] + (bins[j] + 1) * (max_X[pos] - min_X[pos]) / init_bins_count[j];
          start_params[pos] = urand(val_from, val_to);
        }
      }

      Optimizer *opt = new Adam();
      Block adam_settings;
      adam_settings.add_arr("initial_params", start_params);
      adam_settings.add_double("learning_rate", local_search_learning_rate);
      adam_settings.add_int("iterations", local_search_iterations);
      adam_settings.add_bool("verbose", verbose);
      adam_settings.add_arr("derivatives_mult", derivatives_mult);
      opt->optimize(F, min_X, max_X, adam_settings);
      
      std::vector<float> local_best_params;
      float local_best_result;
      local_best_params = opt->get_best_result(&local_best_result);
      if (local_best_result < best_result)
      {
        best_result = local_best_result;
        best_params = local_best_params;
      }
      if (verbose)
      {
        debug("GSA: local search %d. Value: %.4f, best: %.4f\n", b, local_best_result, best_result);
      }

      delete opt;
      b++;
    }
  }
}