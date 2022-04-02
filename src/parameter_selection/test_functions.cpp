#include "generic_optimization_algorithm.h"
#include "common_utils/utility.h"

namespace my_opt
{
void get_test_function(std::string name, my_opt::OptFunction &optF, std::vector<float> &parList_f)
{
    optF.name = name;
    optF.version = 1;
    if (name == "sphere100")
    {
        parList_f = std::vector<float>(100,0);
        optF.f = [&](std::vector<std::vector<float>> &params) -> std::vector<float> {
            std::vector<float> res;
            for (auto &x : params)
            {
                res.push_back(0);
                for (auto &x_i : x)
                    res.back() += MAX(0, 10-x_i*x_i);
            }
            return res;
        };
    }
    else if (name == "cross_in_tray")
    {
        parList_f = std::vector<float>(2,0);
        optF.f = [&](std::vector<std::vector<float>> &params) -> std::vector<float> {
            std::vector<float> res;
            for (auto &x : params)
            {
                x[0] = 20*x[0] - 10;
                x[1] = 20*x[1] - 10;
                res.push_back(0);
                res.back() = 0.0001*pow(abs(sin(x[0])*sin(x[1])*exp(abs(100 - sqrt(x[0]*x[0] + x[1]*x[1])/PI)))+1,0.1);
            }
            return res;
        };
    }
    else if (name == "schaffer_2")
    {
        parList_f = std::vector<float>(2,0);
        optF.f = [&](std::vector<std::vector<float>> &params) -> std::vector<float> {
            std::vector<float> res;
            for (auto &x : params)
            {
                x[0] = 200*x[0] - 100;
                x[1] = 200*x[1] - 100;
                res.push_back(0);
                res.back() = 10 - (0.5 + (SQR(sin(x[0]*x[0] - x[1]*x[1])) - 0.5)/SQR(1 + 0.001*(x[0]*x[0] + x[1]*x[1])));
            }
            return res;
        };
    }
    else if (name == "styblinski_tang10")
    {
        parList_f = std::vector<float>(10,0);
        optF.f = [&](std::vector<std::vector<float>> &params) -> std::vector<float> {
            std::vector<float> res;
            for (auto &x : params)
            {
                res.push_back(0);
                double sum = 0;
                for (auto &x_i : x)
                {
                    x_i = 10*x_i - 5;
                    sum += x_i*x_i*x_i*x_i - 16*x_i*x_i + 5*x_i;
                }
                res.back() = -sum / (2*10);
            }
            return res;
        };
    }
    else if (name == "styblinski_tang100")
    {
        parList_f = std::vector<float>(100,0);
        optF.f = [&](std::vector<std::vector<float>> &params) -> std::vector<float> {
            float max_value = 39.166;
            std::vector<float> res;
            for (auto &x : params)
            {
                res.push_back(0);
                double sum = 0;
                for (auto &x_i : x)
                {
                    x_i = 10*x_i - 5;
                    sum += x_i*x_i*x_i*x_i - 16*x_i*x_i + 5*x_i;
                }
                res.back() = -sum / (2*100*max_value);
            }
            return res;
        };
    }
}
}