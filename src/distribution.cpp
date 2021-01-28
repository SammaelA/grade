#include "distribution.h"
#include <random>
#include <iostream>
#include "tinyEngine/utility.h"
int base_seed = 1234;
Uniform base = Uniform(0,1,base_seed);
UniformInt basei = UniformInt(0,INT_MAX,base_seed);
Normal::Normal(double a, double sigma, int seed)
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
Uniform::Uniform(double from, double to,int seed)
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
UniformInt::UniformInt(double from, double to, int seed)
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
double srand(uint64_t s, uint64_t &x, uint64_t&w, double from, double to)
{  
    x *= x; 
    x += (w += s); 
    x = (x>>32) | (x<<32);
    uint32_t mx = ~0;
    return ((double)(x % mx)/mx)*(to - from) + from;
}
double srandi(uint64_t s, uint64_t &x, uint64_t&w, int from, int to)
{   
    x *= x; 
    x += (w += s); 
    x = (x>>32) | (x<<32);
    return x % (to - from) + from;
}