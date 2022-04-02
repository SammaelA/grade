#pragma once
#include <vector>
#include <string>
#include <functional>

namespace my_opt
{
struct MetaParameters
{
    virtual void RW_vector(std::vector<float> &res) = 0;
};
struct OptFunction
{
    std::string name = "default";
    int version = 0;
    std::function<std::vector<float>(std::vector<std::vector<float>> &)> f;
};
struct ExitConditions
{
    float time_elapsed_seconds = 45*60;
    float function_reached = 1;
    int function_calculated = 1000;
    int generations = 10000;
};
class Optimizer
{
public:
    virtual void perform(std::vector<float> &param_list, MetaParameters *params, ExitConditions exit_conditions,
                 const OptFunction &opt_f,
                 std::vector<std::pair<float,std::vector<float>>> &best_results,
                 std::vector<std::vector<float>> &initial_types) = 0;
    virtual ~Optimizer() {};
};
void get_test_function(std::string name, my_opt::OptFunction &optF, std::vector<float> &parList_f);
};