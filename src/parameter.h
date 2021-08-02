#pragma once

#include "distribution.h"
#include "tinyEngine/utility.h"
#include <vector>
enum RandomnessLevel
{
    NO_RANDOM,
    EXPLICIT_REGENERATION,
    REGENERATE_ON_STATE_CHANGE,
    REGENERATE_ON_GET
};
enum ParameterMaskValues
{
    CONSTANT,
    ONE_VALUE,
    LIST_OF_VALUES,
    FULL
};
enum ParameterVariablesSet
{
    ONLY_BASE_VALUES,
    BASE_VALUES_AND_QS,
    ALL_VALUES
};
struct TreeStructureParameters;
template <typename T>
class Parameter
{
public:
    friend struct TreeStructureParameters;
    float get_min() { return minValue;}
    float get_max() { return maxValue;}
    std::string to_string()
    {
        std::string str = "Parameter ";
        str = str + std::to_string(baseValue);
        str = str + " {";
        for (auto q : state_qs)
        {
            str = str + std::to_string(q*baseValue);
            str = str + " ";
        }
        str = str + "} ";
        str = str + std::to_string(a) + " ";
        str = str + std::to_string(sigma) + " ";
        str = str + std::to_string(from) + " ";
        str = str + std::to_string(to) + " ";
        str = str + std::to_string(normal_part);

        return str;
    }   
    void random_regenerate()
    {
        if (normal_part < 1e-4)
            randomValue = uniform->get();
        else if (1 - normal_part < 1e-4)
            randomValue = normal->get();
        else 
            randomValue = (T)(normal_part*normal->get() + (1-normal_part)*uniform->get());
    }
    Parameter &operator=(const T &base_value)
    {
        this->baseValue = baseValue;
        randomnessLevel = NO_RANDOM;
        return *this;
    }
    Parameter &operator=(const Parameter &par)
    {
        baseValue = par.baseValue;
        randomValue = par.randomValue;
        maxValue = par.maxValue;
        minValue = par.minValue;
        minMaxDefined = par.minMaxDefined;
        state_qs = par.state_qs;
        state = par.state;
        normal = par.normal;
        uniform = par.uniform;
        randomnessLevel = par.randomnessLevel;
        a = par.a;
        sigma = par.sigma;
        from = par.from;
        to = par.to;
        normal_part = par.normal_part;
        return *this;
    }
    void set_no_override_minmax(const Parameter &par)
    {
        T mn,mx;
        bool mmdef = minMaxDefined || par.minMaxDefined;
        if (false && !minMaxDefined && par.minMaxDefined)
        {
            mn = par.minValue;
            mx = par.maxValue;
        }
        else
        {
            mn = minValue;
            mx = maxValue;
        }
        
        *this = par;
        minValue = mn;
        maxValue = mx;
        minMaxDefined = mmdef;
    }
    T get()
    {
        if (randomnessLevel == NO_RANDOM)
        {
            if (state == -1)
                return baseValue;
            else
                return (T)(state_qs[state]*baseValue);
        }
        if (randomnessLevel == REGENERATE_ON_GET)
            random_regenerate();
        T value = baseValue + randomValue;
        if (state != -1)
            value *= state_qs[state];

        if (value > maxValue)
            value = maxValue;
        else if (value < minValue)
            value = minValue;
        return value;
    }
    T operator()()
    {
        return get();
    }
    void set_state(int state)
    {
        if (state < 0 || (uint)state >= state_qs.size())
            state = -1;
        if (randomnessLevel == REGENERATE_ON_STATE_CHANGE && this->state != state)
            random_regenerate();
        this->state = state;
    }
    Parameter(const Parameter &par)
    {
        baseValue = par.baseValue;
        randomValue = par.randomValue;
        maxValue = par.maxValue;
        minValue = par.minValue;
        minMaxDefined = par.minMaxDefined;
        state_qs = par.state_qs;
        state = par.state;
        normal = par.normal;
        uniform = par.uniform;
        randomnessLevel = par.randomnessLevel;
        a = par.a;
        sigma = par.sigma;
        from = par.from;
        to = par.to;
        normal_part = par.normal_part;
    }
    Parameter(T base, T minValue = (T)(-1e10), T maxValue = (T)(1e10))
    {
        if (minValue > maxValue)
        {
            float t = minValue;
            minValue = maxValue;
            maxValue = t;
        }
        baseValue = CLAMP(base, minValue, maxValue);
        randomnessLevel = NO_RANDOM;
        this->maxValue = maxValue;
        this->minValue = minValue;
        minMaxDefined = (minValue > -1e10 && maxValue < 1e10);
        state = -1;
    }
    Parameter(T base, std::vector<T> stateParams,
              T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : Parameter(base, minValue, maxValue)
    {
        if (minValue > maxValue)
        {
            float t = minValue;
            minValue = maxValue;
            maxValue = t;
        }
        for (T &par : stateParams)
        {
            state_qs.push_back(base > 0 ? ((double)CLAMP(par, minValue, maxValue)/base) : 1);
        }
    }
    Parameter(T base, std::vector<T> stateParams, RandomnessLevel rand_level, Normal *randomizer,
              T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : Parameter(base, stateParams, minValue, maxValue)
    {
        this->randomnessLevel = rand_level;
        this->normal = randomizer;
        if (!randomizer)
            rand_level = NO_RANDOM;
        else
        {
            a = randomizer->get_a();
            sigma = randomizer->get_sigma();
        }
        from = 0;
        to = 0;
        normal_part = 1;

        if (rand_level != NO_RANDOM)
            random_regenerate();
    }
    Parameter(T base, RandomnessLevel rand_level, Normal *randomizer,
              T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : Parameter(base, minValue, maxValue)
    {
        this->randomnessLevel = rand_level;
        this->normal = randomizer;
        if (!randomizer)
            rand_level = NO_RANDOM;
        else
        {
            a = randomizer->get_a();
            sigma = randomizer->get_sigma();
        }

        a = randomizer->get_a();
        sigma = randomizer->get_sigma();
        from = 0;
        to = 0;
        normal_part = 1;
        
        if (rand_level != NO_RANDOM)
            random_regenerate();

    }
        Parameter(T base, std::vector<T> stateParams, RandomnessLevel rand_level, Uniform *randomizer,
              T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : Parameter(base, stateParams, minValue, maxValue)
    {
        this->randomnessLevel = rand_level;
        this->uniform = randomizer;
        if (!randomizer)
            rand_level = NO_RANDOM;
        else
        {
            from = randomizer->get_from();
            to = randomizer->get_to();
        }
        a = 0;
        sigma = 0;
        normal_part = 0;
        
        if (rand_level != NO_RANDOM)
            random_regenerate();
    }
    Parameter(T base, RandomnessLevel rand_level, Uniform *randomizer,
              T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : Parameter(base, minValue, maxValue)
    {
        this->randomnessLevel = rand_level;
        this->uniform = randomizer;
        if (!randomizer)
            rand_level = NO_RANDOM;
        else
        {
            from = randomizer->get_from();
            to = randomizer->get_to();
        }      
        a = 0;
        sigma = 0;
        normal_part = 0;

        if (rand_level != NO_RANDOM)
            random_regenerate();
    }
    Parameter(T base, T min_val, T max_val, std::vector<float> state_qs, float a, float sigma, float from, float to,
              float normal_part, RandomnessLevel rand_level, Normal *normal = nullptr, Uniform *uniform = nullptr)
    {
        minValue = min_val;
        maxValue = max_val;
        if (minValue > maxValue)
        {
            float t = minValue;
            minValue = maxValue;
            maxValue = t;
        }
        baseValue = CLAMP(base, minValue, maxValue);
        minMaxDefined = (minValue > -1e10 && maxValue < 1e10);
        this->state_qs = state_qs;
        for (float &q : this->state_qs)
        {
            T val = CLAMP(base*q, minValue, maxValue);
            q = baseValue > 0 ? val/baseValue : 1;
        }
        this->a = a;
        this->sigma = sigma;
        this->from = from;
        this->to = to;
        this->normal_part = normal_part;
        this->randomnessLevel = randomnessLevel;
        this->normal = normal;
        this->uniform = uniform;
        if (randomnessLevel != RandomnessLevel::NO_RANDOM && normal_part > 1e-4 && !normal)
        {
            normal = distibutionGenerator.get_normal(a,sigma);
        }
        if (randomnessLevel != RandomnessLevel::NO_RANDOM && 1 - normal_part > 1e-4 && !uniform)
        {
            uniform = distibutionGenerator.get_uniform(from,to);
        } 
    }
private:
    T baseValue;
    T randomValue = 0;
    T maxValue = 1e10;
    T minValue = -1e10;
    bool minMaxDefined = false;
    std::vector<float> state_qs;
    float a = 0,sigma = 0,from = 0,to = 0,normal_part = 0;
    int state = -1;
    Distribution *normal = nullptr;
    Distribution *uniform = nullptr;
    RandomnessLevel randomnessLevel = NO_RANDOM;
};
struct ParameterDesc
{
    ParameterMaskValues mask;
    RandomnessLevel randomnessLevel;
    int var_count;
    Parameter<float> original;
    std::string name;
    ParameterDesc(Parameter<float> &p):original(p) {};
};
struct ParameterTinyDesc
{
    ParameterMaskValues val;
    std::string name;
};
struct TreeStructureParameters
{
    static Parameter<int> from_float(Parameter<float> source);
    static Parameter<float> from_int(Parameter<int> source);
    Parameter<float> max_depth;
    Parameter<float> max_segments;
    Parameter<float> max_branching;
    Parameter<float> growth_iterations;

    Parameter<float> scale;
    Parameter<float> seg_len_mult;
    Parameter<float> leaf_size_mult;
    Parameter<float> base_r;
    Parameter<float> r_split_save_pow;

    Parameter<float> dir_conserv;
    Parameter<float> plane_conserv;
    Parameter<float> spread;
    Parameter<float> phototrop;
    Parameter<float> gravitrop;
    Parameter<float> dir_random;
    Parameter<float> base_angle;
    Parameter<float> base_angle_q;

    Parameter<float> seg_dir_conserv;
    Parameter<float> seg_plane_conserv;
    Parameter<float> seg_spread;
    Parameter<float> seg_phototrop;
    Parameter<float> seg_gravitrop;
    Parameter<float> seg_dir_random;
    Parameter<float> seg_bend;
    Parameter<float> seg_bend_pow;

    Parameter<float> base_branch_feed;
    Parameter<float> base_seg_feed;
    Parameter<float> feed_distribution_min_weight;
    Parameter<float> feed_distribution_d_weight;
    Parameter<float> top_growth_bonus;

    Parameter<float> light_precision;
    Parameter<float> branch_removal;
    Parameter<float> branch_grow_decrease_q;
    Parameter<float> segment_grow_decrease_q;

    Parameter<float> min_branching_chance;
    Parameter<float> max_branching_chance;
    Parameter<float> branching_power;

    Parameter<float> r_deformation_levels;
    Parameter<float> r_deformation_points;
    Parameter<float> r_deformation_power;
    void set_state(int state)
    {
        max_depth.set_state(state);
        max_segments.set_state(state);
        max_branching.set_state(state);
        growth_iterations.set_state(state);

        scale.set_state(state);
        seg_len_mult.set_state(state);
        leaf_size_mult.set_state(state);
        base_r.set_state(state);
        r_split_save_pow.set_state(state);

        dir_conserv.set_state(state);
        plane_conserv.set_state(state);
        spread.set_state(state);
        phototrop.set_state(state);
        gravitrop.set_state(state);
        dir_random.set_state(state);
        base_angle.set_state(state);
        base_angle_q.set_state(state);

        seg_dir_conserv.set_state(state);
        seg_plane_conserv.set_state(state);
        seg_spread.set_state(state);
        seg_phototrop.set_state(state);
        seg_gravitrop.set_state(state);
        seg_dir_random.set_state(state);
        seg_bend.set_state(state);
        seg_bend_pow.set_state(state);

        base_branch_feed.set_state(state);
        base_seg_feed.set_state(state);
        feed_distribution_min_weight.set_state(state);
        feed_distribution_d_weight.set_state(state);
        top_growth_bonus.set_state(state);

        light_precision.set_state(state);
        branch_removal.set_state(state);
        branch_grow_decrease_q.set_state(state);
        segment_grow_decrease_q.set_state(state);

        min_branching_chance.set_state(state);
        max_branching_chance.set_state(state);
        branching_power.set_state(state);

        r_deformation_levels.set_state(state);
        r_deformation_points.set_state(state);
        r_deformation_power.set_state(state);
    }

    TreeStructureParameters() : max_depth(4,4,4),
                                max_segments(40, std::vector<float>{40, 14, 14, 7, 7}, 10, 40),
                                max_branching(1,1,1),
                                growth_iterations(200,200,200),

                                scale(3,3,3),
                                seg_len_mult(2.25, std::vector<float>{2.25, 1.75, 1, 0.55, 0.4}, 0.1, 5),
                                leaf_size_mult(2.25),
                                base_r(0.9, std::vector<float>{0.9, 0.75, 0.35, 0.15, 0.12}, 0.5, 1.5),
                                r_split_save_pow(3, std::vector<float>{1.5, 1.7, 2.3, 3, 2}, REGENERATE_ON_GET, distibutionGenerator.get_normal(0, 0.25), 1, 4),

                                dir_conserv(3, std::vector<float>{2, 2, 3, 2, 2}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.4, 0.4), 0.1, 10),
                                plane_conserv(3, std::vector<float>{2, 2, 3, 2, 2}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.4, 0.4), 0.1, 10),
                                spread(3, std::vector<float>{1, 3, 2, 2, 2}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.25, 25), 0.1, 10),
                                phototrop(5, std::vector<float>{5, 3, 2, 0, 0}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.2, 0.2), 0.1, 10),
                                gravitrop(5, std::vector<float>{5, 5, 2, 0.75, 0.25}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.1, 0.1), 0.1, 10),
                                dir_random(1, std::vector<float>{1, 1, 1, 1, 1}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.4, 0.4), 0.1, 10),
                                base_angle(3.141/4, 3.141/4, 3.141/4),
                                base_angle_q(3, std::vector<float>{1, 3, 2, 2, 2}, 0.1, 10),

                                seg_dir_conserv(20, std::vector<float>{10, 10, 10, 10, 20}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.05, 0.05), 0, 25),
                                seg_plane_conserv(20, std::vector<float>{10, 10, 10, 10, 20}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.05, 0.05), 0, 25),
                                seg_spread(1, std::vector<float>{1, 1, 0.5, 0.5, 1}, 0.1, 5),
                                seg_phototrop(7, std::vector<float>{7, 3.5, 1.5, 0.4, 0.1}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-1, 2), 0, 50),
                                seg_gravitrop(2, std::vector<float>{0.5, 0.2, 0.07, 0.04, 0.5}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-2, 4), 0, 50),
                                seg_dir_random(1, std::vector<float>{0.5, 1, 2, 3, 5}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-1, 1), 0, 50),
                                seg_bend(1, std::vector<float>{1, 1, 2, 3, 5}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-1, 1), 0, 50),
                                seg_bend_pow(2, std::vector<float>{2, 2, 2, 2}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-0.05, 0.05), 0, 50),

                                base_branch_feed(800, std::vector<float>{800, 350, 200, 40, 40}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-200, 200), 1, 1000),
                                base_seg_feed(200, std::vector<float>{200, 150, 120, 40, 30}, REGENERATE_ON_GET, distibutionGenerator.get_uniform(-30, 30), 1, 250),
                                feed_distribution_min_weight(0.07, 0.07, 0.07),
                                feed_distribution_d_weight(0.05, 0.05, 0.05),
                                top_growth_bonus(0.0, 0, 0),

                                light_precision(0.25, 0.25, 0.25),
                                branch_removal(1.2, std::vector<float>{0, 0.75, 1, 1.2, 1.2}, 0.0, 1.5),
                                branch_grow_decrease_q(0.5, 0.5, 0.5),
                                segment_grow_decrease_q(0.05, 0.05, 0.05),

                                min_branching_chance(0.7, std::vector<float>{0.0, 0.1, 0.5, 0.75, 0.7}, 0, 1),
                                max_branching_chance(1, 1, 1),
                                branching_power(1.2, std::vector<float>{1.2, 0.8, 0.5, 0.5, 0.4}, 0.5, 2.5),
                                r_deformation_levels(2, 2, 2),
                                r_deformation_points(8,std::vector<float>{8,3}, 8, 8),
                                r_deformation_power(0, std::vector<float>{0,0}, REGENERATE_ON_GET, distibutionGenerator.get_normal(0,0.025), 0, 1)
    {
    }
    void get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                            ParameterVariablesSet v_set = ParameterVariablesSet::ALL_VALUES);
    void get_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, 
                           ParameterVariablesSet v_set = ParameterVariablesSet::ONLY_BASE_VALUES);
    void load_from_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data,
                                 ParameterVariablesSet v_set = ParameterVariablesSet::ONLY_BASE_VALUES);
};
class ParameterSetWrapper
{
    std::vector<TreeStructureParameters> params;
    int state = 0;
    public:
    ParameterSetWrapper()
    {
        params.push_back(TreeStructureParameters());
    }
    ParameterSetWrapper(TreeStructureParameters base, int n_params)
    {
        for (int i=0;i<n_params;i++)
        {
            params.push_back(base);
            params.back().set_state(i);
        }
    }
    void set_state(int _state) {state = _state;}
    TreeStructureParameters &operator()() {return params[state];}
};