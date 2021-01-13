#include "distribution.h"
#include <random>
#include <iostream>
#include "tinyEngine/utility.h"
int base_seed = 1234;
Uniform base = Uniform(0,1);
UniformInt basei = UniformInt(0,INT_MAX);
Normal::Normal(double a, double sigma)
{
    this->a = a;
    this->sigma = sigma;
    d = std::normal_distribution<double>(a, sigma);
    gen.seed(base_seed);
}
double Normal::get()
{
    return d(gen);
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
Uniform::Uniform(double from, double to)
{
    this->from = from;
    this->to = to;
    delta = to - from;
    if (delta < 0)
        logerr("wrong interval for uniform distribution [%lg %lg)", from, to);
    d = std::uniform_real_distribution<double>(from,to);
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
UniformInt::UniformInt(double from, double to)
{
    this->from = from;
    this->to = to;
    delta = to - from;
    if (delta < 0)
        logerr("wrong interval for uniform distribution [%lg %lg)", from, to);
    d = std::uniform_int_distribution<long>(from,to);
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
    return base.get()*(to - from) + from;
}
double urandi(int from, int to)
{
    return basei.geti()%(to - from) + from;
}