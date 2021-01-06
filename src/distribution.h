#pragma once
#include <random>
class Distribution
{
public:
    Distribution(){};
    virtual double get() = 0;
    virtual double *get_series(unsigned size) = 0;
};
class Normal : public Distribution
{
    double a;
    double sigma;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d;

public:
    Normal(double a = 0, double sigma = 1);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};
class Uniform : public Distribution
{
    double from;
    double to;
    double delta;

public:
    Uniform(double from = 0, double to = 1);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};