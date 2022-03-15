#pragma once

#include "common_utils/distribution.h"
#include "common_utils/utility.h"
#include "save_utils/blk.h"

#include <vector>
DEFINE_ENUM_WITH_STRING_CONVERSIONS(RandomnessLevel,(NO_RANDOM)(EXPLICIT_REGENERATION)(REGENERATE_ON_STATE_CHANGE)(REGENERATE_ON_GET))

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
struct ParametersSet;
struct WeberPennParameters;
struct ParameterList;
template <typename T>
class Parameter
{
public:
    friend struct ParametersSet;
    friend struct WeberPennParameters;
    
    float get_base() { return baseValue;} 
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
    explicit Parameter(std::vector<T> stateParams,
                       T minValue = (T)(-1e10), T maxValue = (T)(1e10)) : 
                       Parameter(stateParams.front(),stateParams, minValue, maxValue)
    {

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
protected:
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
struct ParametersSet
{
    static Parameter<int> from_float(Parameter<float> source);
    static Parameter<float> from_int(Parameter<int> source);
    virtual void get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                            ParameterVariablesSet v_set = ParameterVariablesSet::ALL_VALUES) {};
    virtual void get_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, 
                           ParameterVariablesSet v_set = ParameterVariablesSet::ONLY_BASE_VALUES);
    virtual void load_from_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data,
                                 ParameterVariablesSet v_set = ParameterVariablesSet::ONLY_BASE_VALUES);
    virtual void save_to_blk(Block &b);
    virtual void load_from_blk(Block &b);
    virtual void set_state(int state) {};
    virtual glm::vec3 get_tree_max_size() = 0;
    virtual ParametersSet *copy() { return nullptr;};
    virtual float get_scale_factor() {return 1;}
    virtual void write_parameter_list(ParameterList &list) { RW_parameter_list(true, list);}
    virtual void read_parameter_list(ParameterList &list) {RW_parameter_list(false, list);}
private:
    virtual void RW_parameter_list(bool write, ParameterList &list){};
};

struct CategorialParameter
{
    int val = 0;
    std::vector<int> possible_values;
    operator int () {return val;}
    CategorialParameter(int _val) 
    {
        val = _val;
        possible_values = {_val};
    }
    bool fixed() { return possible_values.size() == 1;}
};

struct OrdinalParameter
{
    int val = 0;
    int min_val = 0;
    int max_val = 0;
    operator int () {return val;}
    OrdinalParameter(int _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    bool fixed() { return max_val <= min_val;}
};

struct ContinuousParameter
{
    float val = 0;
    float min_val = 0;
    float max_val = 0;
    operator float () {return val;}
    ContinuousParameter(int _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    ContinuousParameter(float _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    bool fixed() { return max_val <= min_val;}
};

struct ParameterList
{
    std::map<std::string, CategorialParameter> categorialParameters;
    std::map<std::string, OrdinalParameter> ordinalParameters;
    std::map<std::string, ContinuousParameter> continuousParameters;

    void print();
    void load_borders_from_blk(Block &b);
    void to_simple_list(std::vector<float> &list, bool normalized = false, bool remove_fixed_params = false);
    void from_simple_list(std::vector<float> &list, bool normalized = false, bool remove_fixed_params = false);
    float diff(ParameterList &list, bool normalized = false, bool remove_fixed_params = false);
};