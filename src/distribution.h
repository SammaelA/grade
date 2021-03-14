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
    Normal(double a = 0, double sigma = 1, int seed = 12345);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};
class Uniform : public Distribution
{
    double from;
    double to;
    double delta;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> d;

public:
    Uniform(double from = 0, double to = 1, int seed = 12345);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};
class UniformInt : public Distribution
{
    double from;
    double to;
    double delta;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<long> d;

public:
    UniformInt(double from = 0, double to = 1, int seed = 12345);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
    long geti();
    long *get_seriesi(unsigned size);
};
class DiscreteGeneral : public Distribution
{
public:
    DiscreteGeneral();
    DiscreteGeneral(const std::vector<double> &values, const std::vector<double> &weights, int seed = 12345);
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
private:
    Uniform u;
    std::vector<double> cumulative_prob;
    std::vector<double> values;
};
double urand(double from = 0.0, double to = 1.0);
double urandi(int from = 0, int to = 1);
double srand(uint64_t seed, uint64_t &x, uint64_t&w, double from = 0.0, double to = 1.0);
double srandi(uint64_t seed, uint64_t &x, uint64_t&w, int from = 0, int to = 1);