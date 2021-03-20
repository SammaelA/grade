#include "distribution.h"
#include <random>
#include <iostream>
#include "tinyEngine/utility.h"
int base_seed = 1234;
Uniform base = Uniform(0, 1, base_seed);
UniformInt basei = UniformInt(0, INT_MAX, base_seed);
Normal::Normal(double a, double sigma, int seed)
{
    this->a = a;
    this->sigma = sigma;
    d = std::uniform_real_distribution<double>(0,1);
    gen.seed(base_seed);
}
double Normal::get()
{
    //Box–Muller transform
    double U1 = d(gen);
    double U2 = d(gen);
    
    return sqrt(-2*log(U1))*cos(2*PI*U2)*sigma + a;
}
double *Normal::get_series(unsigned size)
{
    double *p = new double[size];
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = d(gen);
    }
    return p;
}
Uniform::Uniform(double from, double to, int seed)
{
    this->from = from;
    this->to = to;
    delta = to - from;
    if (delta < 0)
        logerr("wrong interval for uniform distribution [%lg %lg)", from, to);
    d = std::uniform_real_distribution<double>(from, to);
    gen.seed(base_seed);
}
double Uniform::get()
{
    return d(gen);
}
double *Uniform::get_series(unsigned size)
{
    double *p = new double[size];
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = get();
    }
    return p;
}
UniformInt::UniformInt(double from, double to, int seed)
{
    this->from = from;
    this->to = to;
    delta = to - from;
    if (delta < 0)
        logerr("wrong interval for uniform distribution [%lg %lg)", from, to);
    d = std::uniform_int_distribution<long>(from, to);
    gen.seed(base_seed);
}
long UniformInt::geti()
{
    return d(gen);
}
long *UniformInt::get_seriesi(unsigned size)
{
    long *p = new long[size];
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = get();
    }
    return p;
}
double UniformInt::get()
{
    return d(gen);
}
double *UniformInt::get_series(unsigned size)
{
    double *p = new double[size];
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = get();
    }
    return p;
}
double urand(double from, double to)
{
    return base.get() * (to - from) + from;
}
double urandi(int from, int to)
{
    return basei.geti() % (to - from) + from;
}
double srand(uint64_t s, uint64_t &x, uint64_t &w, double from, double to)
{
    x *= x;
    x += (w += s);
    x = (x >> 32) | (x << 32);
    uint32_t mx = ~0;
    return ((double)(x % mx) / mx) * (to - from) + from;
}
double srandi(uint64_t s, uint64_t &x, uint64_t &w, int from, int to)
{
    x *= x;
    x += (w += s);
    x = (x >> 32) | (x << 32);
    return x % (to - from) + from;
}
DiscreteGeneral::DiscreteGeneral()
{
    cumulative_prob.push_back(1);
    values.push_back(0);
}
DiscreteGeneral::DiscreteGeneral(std::vector<double> &values, std::vector<double> &weights,
                                 int seed) : u(0, 1, seed)
{
    if (values.size() != weights.size() || values.empty())
    {
        logerr("wrong parameters for general discrete distribution - vectors of values and weights should have same size bigger than zero");
        logerr("%d %d",values.size(),weights.size());
        this->values = {0};
        cumulative_prob.push_back(1);
    }
    else
    {
        this->values = values;
        for (int i = 0; i < weights.size(); i++)
        {
            double w = weights[i] + (i == 0 ? 0 : cumulative_prob.back());
            cumulative_prob.push_back(w);
        }
        int sum = cumulative_prob.back();
        for (int i = 0; i < cumulative_prob.size(); i++)
        {
            cumulative_prob[i] /= sum;
        }
    }
}
double DiscreteGeneral::get()
{
    double rnd = u.get();
    for (int i = 0; i < cumulative_prob.size(); i++)
    {
        if (rnd < cumulative_prob[i])
            return values[i];
        else
            rnd -= cumulative_prob[i];
    }
    return values.back();
}
double *DiscreteGeneral::get_series(unsigned size)
{
}