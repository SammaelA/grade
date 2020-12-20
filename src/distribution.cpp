#include "distribution.h"
#include <random>
#include <iostream>

    Normal::Normal(double a, double sigma)
    {
        this->a = a;
        this->sigma = sigma;
        d = std::normal_distribution<double>(a,sigma);
    }
    double Normal::get()
    {
        return d(gen);
    }
    double *Normal::get_series(unsigned size)
    {
        double *p = new double[size];
        for (unsigned i =0;i<size;i++)
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
        if (delta<0)
            fprintf(stderr,"wrong interval for uniform distribution [%lg %lg)",from,to);
    }
    double Uniform::get()
    {
        return (double)rand()*delta/RAND_MAX + from;
    }
    double *Uniform::get_series(unsigned size)
    {
        double *p = new double[size];
        for (unsigned i =0;i<size;i++)
        {
            p[i] = get();
        }
        return p;
    }