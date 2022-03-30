#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include "generic_optimization_algorithm.h"
#include <chrono>

class SimulatedAnnealing : public my_opt::Optimizer
{
public:
    struct MetaParameters : public my_opt::MetaParameters
    {
        float t_start = 1;
        float t_pow = 1;
        float tries = 1;
        virtual void RW_vector(std::vector<float> &res) {};
    };

    void perform(std::vector<float> &param_list, my_opt::MetaParameters *params, my_opt::ExitConditions exit_conditions,
                 const my_opt::OptFunction &my_opt_f,
                 std::vector<std::pair<float,std::vector<float>>> &best_results,
                 std::vector<std::vector<float>> &initial_types);
private:
    my_opt::OptFunction function;
    my_opt::ExitConditions exitConditions;
    MetaParameters metaParams;
    int free_parameters_cnt = 0;
    std::chrono::steady_clock::time_point t_start, iter_start;
    float best_metric_ever = 0;
    float time_spent_sec = 0; 
    int func_called = 0;
    int iteration_n = 0;

    typedef std::vector<float> State;
    float temp(int iter, int max_iter);
    float F(State &state);
    State next_state(State &prev_state, float mutation_power);
    State initial_state();
    bool should_exit();
};