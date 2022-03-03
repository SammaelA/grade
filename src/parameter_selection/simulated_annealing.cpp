#include "simulated_annealing.h"

void SimulatedAnnealing::perform(ParameterList &param_list, MetaParameters params, ExitConditions exit_conditions,
                                 const std::function<std::vector<float>(std::vector<ParameterList> &)> &f,
                                 std::vector<std::pair<float, ParameterList>> &best_results,
                                 std::vector<ParameterList> &initial_types)
{
    t_start = std::chrono::steady_clock::now();

    function = f;
    original_param_list = param_list;
    metaParams = params;
    exitConditions = exit_conditions;
    free_parameters_cnt = 0;

    for (auto &p : original_param_list.categorialParameters)
    {
        parametersMask.categorialParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }
    for (auto &p : original_param_list.ordinalParameters)
    {
        parametersMask.ordinalParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }
    for (auto &p : original_param_list.continuousParameters)
    {
        parametersMask.continuousParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }

    int iters = exitConditions.function_calculated/metaParams.tries;
    State best_state_ever;
    float best_metric_ever = -1;
    for (int tr = 0;tr<metaParams.tries;tr++)
    {
        int iter = 0;
        State st = initial_state();
        float st_metr = F(st);
        State best_state = st;
        float best_metric = st_metr;
        while (!should_exit() && iter < iters)
        {
            State new_st = next_state(st, 1 - (0.5*(float)(iter)/iters));
            float new_metr = F(new_st);
            if (new_metr > st_metr)
            {
                st_metr = new_metr;
                st = new_st;
            }
            else
            {
                float rnd = urand();
                float diff = st_metr - new_metr;
                float T = temp(iter, iters);
                if (rnd < exp(-diff/T))
                {
                    st_metr = new_metr;
                    st = new_st;                    
                }
                else
                {

                }
            }
            if (st_metr > best_metric)
            {
                best_metric = st_metr;
                best_state = st;
            }
            if (iter % 50 == 0)
                logerr("[%d][%d] val = %.4f best = %.4f", tr,iter, new_metr,best_metric);
            iter++;
        }

        if (best_metric > best_metric_ever)
        {
            best_metric_ever = best_metric;
            best_state_ever = best_state;
        }

        best_results.emplace_back();
        best_results.back().first = best_metric;
        best_results.back().second = original_param_list;
        best_results.back().second.from_simple_list(best_state);
    }
}


float SimulatedAnnealing::temp(int iter, int max_iter)
{
    return metaParams.t_start*pow(1 - (float)iter/max_iter,metaParams.t_pow);
}

float SimulatedAnnealing::F(State &state)
{
    std::vector<ParameterList> type = {original_param_list};
    type.back().from_simple_list(state);
    auto res = function(type);
    return res[0];
}

SimulatedAnnealing::State SimulatedAnnealing::next_state(State &prev_state, float mutation_power)
{
    int mutation_genes_count = 7;
    State G = prev_state;

    for (int gene = 0; gene < mutation_genes_count; gene++)
    {
        bool found = false;
        int tries = 0;
        while (!found && tries < 100)
        {
            int pos = urandi(0, G.size());
            int g_pos = pos;
            if (pos < parametersMask.categorialParameters.size())
            {
                if (!parametersMask.categorialParameters[pos].second.fixed())
                {
                    if (urand() < mutation_power)
                    {
                        int val_pos = urandi(0, parametersMask.categorialParameters[pos].second.possible_values.size());
                        G[g_pos] = parametersMask.categorialParameters[pos].second.possible_values[val_pos];
                    }
                    found = true;
                }
            }
            else if (pos < parametersMask.categorialParameters.size() + parametersMask.ordinalParameters.size())
            {
                pos -= parametersMask.categorialParameters.size();
                if (!parametersMask.ordinalParameters[pos].second.fixed())
                {
                    float len = mutation_power*(parametersMask.ordinalParameters[pos].second.max_val - parametersMask.ordinalParameters[pos].second.min_val);
                    float from = CLAMP(G[g_pos] - len, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    float to = CLAMP(G[g_pos] + len, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    float val = urand(from,to);
                    val = CLAMP(val, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    G[g_pos] = val;
                    found = true;
                }
            }
            else
            {
                pos -= parametersMask.categorialParameters.size() + parametersMask.ordinalParameters.size();
                //logerr("continuous pos %d mutation", pos);
                if (!parametersMask.continuousParameters[pos].second.fixed())
                {
                    float len = mutation_power*(parametersMask.continuousParameters[pos].second.max_val - parametersMask.continuousParameters[pos].second.min_val);
                    float from = CLAMP(G[g_pos] - len, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    float to = CLAMP(G[g_pos] + len, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    float val = urand(from,to);
                    val = CLAMP(val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    G[g_pos] = val;
                    //logerr("applied pos val %d %f in [%f %f]",pos, val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    found = true;
                }
            }
            tries++;
        }
    }
    return G;
}

SimulatedAnnealing::State SimulatedAnnealing::initial_state()
{
    State g;
    for (auto &p : parametersMask.categorialParameters)
        g.push_back(p.second.possible_values[(int)urandi(0, p.second.possible_values.size())]);
    for (auto &p : parametersMask.ordinalParameters)
        g.push_back((int)(p.second.min_val + urand()*(p.second.max_val - p.second.min_val)));
    for (auto &p : parametersMask.continuousParameters)
        g.push_back((p.second.min_val + urand()*(p.second.max_val - p.second.min_val)));
    return g;
}

bool SimulatedAnnealing::should_exit()
{
    iter_start = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(iter_start - t_start).count();
    time_spent_sec = 0.001*time;
    return (time_spent_sec > exitConditions.time_elapsed_seconds || iteration_n >= exitConditions.generations
    || func_called >= exitConditions.function_calculated || best_metric_ever >= exitConditions.function_reached);
}
