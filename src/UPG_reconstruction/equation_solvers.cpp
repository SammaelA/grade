#include "equation_solvers.h"

void
solver::find_interval(const std::vector<float>& coefs, const float& x1, const float& x2, const int& n, float intervals[2])
{
    int nb = n, nroot = 0;
    float dx = (x2 - x1) / (n * 1000);
    float x = x1;
    float fp = f(coefs, x1, n);

    bool go_further = true;

    for (int i = 0; i < n * 1000 && go_further; i++)
    {
        x += dx;
        float fc = f(coefs, x, n);

        if (fc * fp <= 0.f)
        {
            // std::cout << x - dx << " " << x << std::endl;

            intervals[0] = x - dx;
            intervals[1] = x;

            go_further = false;
        }

        fp = fc;
    }
}

float 
solver::f(const std::vector<float>& coefs, const float& x, const int& n)
{
    float res = 0;

    for (int i = 0; i <= n; i++)
    {
        res += std::pow(x, i) * coefs[i];
    }

    return res;
}

