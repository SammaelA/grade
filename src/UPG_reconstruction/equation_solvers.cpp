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

void
solver::fd(const std::vector<float>& coefs, const float& x, const float& last_x, const int& n, float res[2])
{
    // float delta = 0.0000001;
    // float f1 = 0;

    // for (int i = 0; i <= n; i++)
    // {
    //     res += std::pow(x, i) * coefs[i];
    // }
}

void 
solver::polinomMiltiplier(const float* p1, const float* p2, float* res, const int& n1, const int& n2, const float& polinom_coef)
{
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            res[i + j] += polinom_coef * p1[i] * p2[j];
        }
    }
}

void
solver::coefsDecrease(const std::vector<float>& coefs, const LiteMath::float3& P, const LiteMath::float3& D, float* res)
{
    float res_coefs[10] = {0};

    float x_polinoms[4][4] = {{1, 0, 0, 0}, {P.x, D.x, 0, 0}, {P.x * P.x, 2 * P.x * D.x, D.x * D.x, 0}, {P.x * P.x * P.x, 3 * P.x * P.x * D.x, 3 * P.x * D.x * D.x, D.x * D.x * D.x}};
    float y_polinoms[4][4] = {{1, 0, 0, 0}, {P.y, D.y, 0, 0}, {P.y * P.y, 2 * P.y * D.y, D.y * D.y, 0}, {P.y * P.y * P.y, 3 * P.y * P.y * D.y, 3 * P.y * D.y * D.y, D.y * D.y * D.y}};
    float z_polinoms[4][4] = {{1, 0, 0, 0}, {P.z, D.z, 0, 0}, {P.z * P.z, 2 * P.z * D.z, D.z * D.z, 0}, {P.z * P.z * P.z, 3 * P.z * P.z * D.z, 3 * P.z * D.z * D.z, D.z * D.z * D.z}};

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                float tmp[10] = {0};

                solver::polinomMiltiplier(x_polinoms[i], y_polinoms[j], tmp, i + 1, j + 1, coefs[i + 4 * j + 16 * k]);
                solver::polinomMiltiplier(tmp, z_polinoms[k], res_coefs, i + j + 1, k + 1, 1);
            }
        }
    }
}