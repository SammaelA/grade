#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include <chrono>

class SimulatedAnnealing
{
public:
    struct MetaParameters
    {
        float t_start = 1;
        float t_pow = 1;
        float tries = 1;
    };
    struct ExitConditions
    {
        float time_elapsed_seconds = 45*60;
        float function_reached = 1;
        int function_calculated = 4000;
        int generations = 10000;
    };
    void perform(ParameterList &param_list, MetaParameters params, ExitConditions exit_conditions,
                 const std::function<std::vector<float>(std::vector<ParameterList> &)> &f,
                 std::vector<std::pair<float,ParameterList>> &best_results,
                 std::vector<ParameterList> &initial_types);
private:
    std::function<std::vector<float>(std::vector<ParameterList> &)> function;
    ParameterList original_param_list;
    struct ParametersMask
    {
        std::vector<std::pair<std::string, CategorialParameter>> categorialParameters;
        std::vector<std::pair<std::string, OrdinalParameter>> ordinalParameters;
        std::vector<std::pair<std::string, ContinuousParameter>> continuousParameters;
    } parametersMask;
    ExitConditions exitConditions;
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