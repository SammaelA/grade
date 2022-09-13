#include "common_utils/distribution.h"
#include <random>
#include <iostream>
#include "common_utils/utility.h"
int base_seed = 1234;
DistributionGenerator distibutionGenerator;
#define RND ((float)rand()/RAND_MAX)
Uniform *base = distibutionGenerator.get_uniform(0, 1, base_seed);
UniformInt *basei = distibutionGenerator.get_uniform_int(0, INT_MAX, base_seed);
Normal::Normal(double a, double sigma, int seed)
{
    this->a = a;
    this->sigma = sigma;
}
double Normal::get()
{
    
    //Box–Muller transform
    double U1 = RND;
    double U2 = RND;
    
    return sqrt(-2*log(U1))*cos(2*PI*U2)*sigma + a;
}
double *Normal::get_series(unsigned size)
{
    double *p = safe_new<double>(size, "normal_series");
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = RND;
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
}
double Uniform::get()
{
    return (RND)*(to - from) + from;
}
double *Uniform::get_series(unsigned size)
{
    double *p = safe_new<double>(size, "uniform_series");
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
}
long UniformInt::geti()
{
    return (long)(RND + from) /(to - from);
}
long *UniformInt::get_seriesi(unsigned size)
{
    long *p = safe_new<long>(size, "uniform_int_series");
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = get();
    }
    return p;
}
double UniformInt::get()
{
    return (long)(RND + from) /(to - from);
}
double *UniformInt::get_series(unsigned size)
{
    double *p = safe_new<double>(size, "uniform_int_series_double");
    for (unsigned i = 0; i < size; i++)
    {
        p[i] = get();
    }
    return p;
}
double urand(double from, double to)
{
    return ((float)rand() / RAND_MAX) * (to - from) + from;
}
double urandi(int from, int to)
{
    return to > from ? rand() % (to - from) + from : to;
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
  return nullptr;
}
glm::vec3 rand_dir()
{
    glm::vec3 dir;
    dir.x = urand(-1,1);
    dir.y = urand(-1,1);
    dir.z = urand(-1,1);
    return normalize(dir);
}