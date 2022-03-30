#include "simulated_annealing.h"

void SimulatedAnnealing::perform(std::vector<float> &param_list, my_opt::MetaParameters *params, my_opt::ExitConditions exit_conditions,
                                 const my_opt::OptFunction &my_opt_f,
                                 std::vector<std::pair<float,std::vector<float>>> &best_results,
                                 std::vector<std::vector<float>> &initial_types)
{
    t_start = std::chrono::steady_clock::now();

    function = my_opt_f;
    MetaParameters *ga_p = dynamic_cast<MetaParameters*>(params);
    if (ga_p)
        metaParams = *ga_p;
    exitConditions = exit_conditions;
    free_parameters_cnt = param_list.size();

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
        best_results.back().second = best_state;
    }
}


float SimulatedAnnealing::temp(int iter, int max_iter)
{
    return metaParams.t_start*pow(1 - (float)iter/max_iter,metaParams.t_pow);
}

float SimulatedAnnealing::F(State &state)
{
    std::vector<std::vector<float>> type = {state};
    auto res = function.f(type);
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
            found = true;
            G[g_pos] = CLAMP(G[g_pos] + mutation_power*urand(-1,1),0,1);
            tries++;
        }
    }
    return G;
}

SimulatedAnnealing::State SimulatedAnnealing::initial_state()
{
    State g(free_parameters_cnt, 0);
    for(int i=0;i<free_parameters_cnt;i++)
        g[i] = urand();
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
