#pragma once
#include "generic_optimization_algorithm.h"
#include "common_utils/LiteMath_ext.h"
#include "common_utils/blk.h"

namespace my_opt
{
    struct FunctionStat
    {
        FunctionStat(int variables_cnt = 0);
        void print();
        void save_load_blk(Block &b, bool save);
        static constexpr int Q_NUM = 10;
        std::string name;
        int version;
        int quantiles[Q_NUM+1];
        int tries = 0;
        std::vector<int[Q_NUM][Q_NUM+1]> by_variable_values;
        std::vector<float[Q_NUM]> marks;
        int marks_q[2*Q_NUM][2];
    };
    void get_function_stat(std::vector<float> &param_list, const OptFunction &my_opt_f, FunctionStat &stat, int count, bool reload = false);
};