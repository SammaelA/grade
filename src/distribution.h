#pragma once
#include <random>

class DistributionGenerator;
class DiscreteGeneral;
class Distribution
{
public:
    Distribution(){};
    virtual double get() = 0;
    virtual double *get_series(unsigned size) = 0;
};
class Normal : public Distribution
{
    friend class DistributionGenerator;
    double a;
    double sigma;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> d;
public:
    Normal(double a = 0, double sigma = 1, int seed = 12345);
    double get_a() {return a;}
    double get_sigma() {return sigma;}
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};
class Uniform : public Distribution
{
    friend class DistributionGenerator;
    friend class DiscreteGeneral;

    double from;
    double to;
    double delta;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> d;
public:
    Uniform(double from = 0, double to = 1, int seed = 12345);
    double get_from() {return from;}
    double get_to() {return to;}
    virtual double get() override;
    virtual double *get_series(unsigned size) override;
};
class UniformInt : public Distribution
{
    friend class DistributionGenerator;

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
    friend class DistributionGenerator;
public:
    DiscreteGeneral();
    DiscreteGeneral(std::vector<double> &values, std::vector<double> &weights, int seed = 12345);
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

class DistributionGenerator
{
    std::vector<Normal *> n_allocs;
    std::vector<Uniform *> u_allocs;
    std::vector<UniformInt *> ui_allocs;
    std::vector<DiscreteGeneral *> dg_allocs;
public:
    Normal *get_normal(double a = 0, double sigma = 1, int seed = 12345)
    {
        n_allocs.push_back(new Normal(a,sigma,seed));
        return (Normal *)(n_allocs.back());
    }
    Uniform *get_uniform(double a = 0, double sigma = 1, int seed = 12345)
    {
        u_allocs.push_back(new Uniform(a,sigma,seed));
        return (Uniform *)(u_allocs.back());
    }
    UniformInt *get_uniform_int(double a = 0, double sigma = 1, int seed = 12345)
    {
        ui_allocs.push_back(new UniformInt(a,sigma,seed));
        return (UniformInt *)(ui_allocs.back());
    }
    DiscreteGeneral *get_discrete_general(std::vector<double> &values, std::vector<double> &weights, int seed = 12345)
    {
        dg_allocs.push_back(new DiscreteGeneral(values,weights,seed));
        return (DiscreteGeneral *)(dg_allocs.back());
    }
    void clear()
    {
        #define DEL(allocs)\
        for (int i=0;i<allocs.size();i++)\
        {\
            if (allocs[i])\
            {\
                delete allocs[i];\
                allocs[i] = nullptr;\
            }\
        }\
        allocs.clear();

        DEL(n_allocs);
        DEL(u_allocs);
        DEL(ui_allocs);
        DEL(dg_allocs);
    }
};
extern DistributionGenerator distibutionGenerator;