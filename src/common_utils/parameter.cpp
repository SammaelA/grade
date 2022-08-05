#include "common_utils/parameter.h"
#include "common_utils/utility.h"


Parameter<int> ParametersSet::from_float(Parameter<float> source)
{
    std::vector<float> qs = source.state_qs;
    
    Parameter<int> par((int)source.baseValue,(int)source.minValue,(int)source.maxValue, qs, source.a, source.sigma, 
                       source.from, source.to, source.normal_part, source.randomnessLevel,
                       dynamic_cast<Normal*>(source.normal), dynamic_cast<Uniform*>(source.uniform));
    return par;
}
Parameter<float> ParametersSet::from_int(Parameter<int> source)
{
    std::vector<float> qs = source.state_qs;
    
    Parameter<float> par((float)source.baseValue,(float)source.minValue,(float)source.maxValue, qs, source.a, source.sigma, 
                       source.from, source.to, source.normal_part, source.randomnessLevel,
                       dynamic_cast<Normal*>(source.normal), dynamic_cast<Uniform*>(source.uniform));
    return par;
}

void ParametersSet::save_to_blk(Block &b)
{
    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    this->get_parameter_list(list,ParameterVariablesSet::ALL_VALUES);

    for (auto &p : list)
    {
        Block *par = new Block();
        par->add_double("baseValue",p.second.baseValue);
        par->add_arr("qs",p.second.state_qs);
        par->add_double("minValue",p.second.minValue);
        par->add_double("maxValue",p.second.maxValue);
        par->add_double("a",p.second.a);
        par->add_double("sigma",p.second.sigma);
        par->add_double("from",p.second.from);
        par->add_double("to",p.second.to);
        par->add_double("normal_part",p.second.normal_part);
        par->add_string("randomnessLevel",ToString(p.second.randomnessLevel));
        b.add_block(p.first.name,par);
        delete par;
    }
}

void ParametersSet::load_from_blk(Block &b)
{
    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    this->get_parameter_list(list,ParameterVariablesSet::ALL_VALUES);

    for (auto &p : list)
    {
        Block *par = b.get_block(p.first.name);
        if (!par)
            continue;
        
        float val = par->get_double("baseValue",p.second.baseValue);
        std::vector<float> state_qs;
        par->get_arr("qs",state_qs);
        float minValue = par->get_double("minValue",p.second.minValue);
        float maxValue = par->get_double("maxValue",p.second.maxValue);
        float a = par->get_double("a",p.second.a);
        float sigma = par->get_double("sigma",p.second.sigma);
        float from = par->get_double("from",p.second.from);
        float to = par->get_double("to",p.second.to);
        float normal_part = par->get_double("normal_part",p.second.normal_part);
        std::string rl_str = par->get_string("randomnessLevel",ToString(p.second.randomnessLevel));
        RandomnessLevel rl = RandomnessLevel::NO_RANDOM;
        for (int i=0;i<(int)(RandomnessLevel::REGENERATE_ON_GET);i++)
        {
            if (ToString((RandomnessLevel)i) == rl_str) 
            {
                rl = (RandomnessLevel)i;
                break;
            }
        }
        p.second = Parameter<float>(val,minValue,maxValue,state_qs,a,sigma,from,to,
                                    normal_part,rl,nullptr,nullptr);
    }
}

void ParameterList::print()
{
    debug("====Parameter List===\n");
    debug("Categorials:\n");
    int n = 0;
    int fn = 0;
    for (auto &p : categorialParameters)
    {
        if(p.second.fixed())
        {
            debug("[%d][%d]%s = %d\n",n,fn, p.first.c_str(), p.second.val);
        }
        else
        {
            debug("[%d][%d]%s = %d in {",n,fn, p.first.c_str(), p.second.val);  
            for (auto &v : p.second.possible_values)
                debug("%d ", v);
            debug("}\n");         
            fn++;
        }
        n++;
    }
    debug("Ordinals:\n");
    for (auto &p : ordinalParameters)
    {
        if(p.second.fixed())
        {
            debug("[%d][%d]%s = %d\n", n, fn, p.first.c_str(), p.second.val);
        }
        else
        {
            debug("[%d][%d]%s = %d in [%d %d]\n", n, fn, p.first.c_str(), p.second.val, p.second.min_val, p.second.max_val);    
            fn++;
        }
        n++;    
    }
    debug("Continuous:\n");
    for (auto &p : continuousParameters)
    {
        if(p.second.fixed())
        {
            debug("[%d][%d]%s = %.2f\n",n,fn, p.first.c_str(), p.second.val);
        }
        else
        {
            debug("[%d][%d]%s = %.2f in [%.2f %.2f]\n",n,fn, p.first.c_str(), p.second.val, p.second.min_val, p.second.max_val); 
            fn++;
        }
        n++;       
    }
}

void ParameterList::load_borders_from_blk(Block &b)
{
    for (int i=0;i<b.size();i++)
    {
        if (b.get_type(i) == Block::ValueType::VEC2)
        {
            //min_max for ordinal or continuous
            auto it = ordinalParameters.find(b.get_name(i));
            if (it != ordinalParameters.end())
            {
                glm::vec2 min_max = b.get_vec2(i, glm::vec2(it->second.val, it->second.val));
                it->second.min_val = min_max.x;
                it->second.max_val = min_max.y;
            }
            else
            {
                auto it = continuousParameters.find(b.get_name(i));
                if (it != continuousParameters.end())
                {
                    glm::vec2 min_max = b.get_vec2(i, glm::vec2(it->second.val, it->second.val));
                    it->second.min_val = min_max.x;
                    it->second.max_val = min_max.y;
                }
                else
                {
                    logerr("%s has wrong parameter name or type in borders list", b.get_name(i).c_str());
                }
            }
        }
        else if (b.get_type(i) == Block::ValueType::ARRAY)
        {
            //categorial
            auto it = categorialParameters.find(b.get_name(i));
            if (it != categorialParameters.end())
            {
                b.get_arr(i, it->second.possible_values, true);
            }
            else
            {
                logerr("%s has wrong parameter name or type in borders list", b.get_name(i));
            }
        }
    }
}

void ParameterList::to_simple_list(std::vector<float> &list, bool normalized, bool remove_fixed_params)
{
    list = {};
    if (normalized)
    {
        for (auto &p : categorialParameters)
        {
            if (p.second.fixed())
            {
                if (!remove_fixed_params)
                    list.push_back(0);
            }
            else
            {
                int pos = 0;
                for (int i=0;i<p.second.possible_values.size();i++)
                {
                    if (p.second.possible_values[i] = p.second.val)
                    {
                        pos = i;
                        break;
                    }
                }
                list.push_back((float)pos/p.second.possible_values.size());
            }
        }
        for (auto &p : ordinalParameters)
        {
            if (!p.second.fixed() || !remove_fixed_params)
                list.push_back(p.second.fixed() ? 0 : (float)(p.second.val - p.second.min_val) / (p.second.max_val - p.second.min_val));
        }
        for (auto &p : continuousParameters)
        {
            if (!p.second.fixed() || !remove_fixed_params)
                list.push_back(p.second.fixed() ? 0 : (float)(p.second.val - p.second.min_val) / (p.second.max_val - p.second.min_val));
        }
    }
    else
    {
        for (auto &p : categorialParameters)
            list.push_back(p.second.val);
        for (auto &p : ordinalParameters)
            list.push_back(p.second.val);
        for (auto &p : continuousParameters)
            list.push_back(p.second.val);
    }
}
void ParameterList::from_simple_list(std::vector<float> &list, bool normalized, bool remove_fixed_params)
{
    int n = 0;
    if (normalized)
    {
        for (auto &p : categorialParameters)
        {
            if (!p.second.fixed())
                p.second.val = p.second.possible_values[(int)list[n]*p.second.possible_values.size()];
            if (!p.second.fixed() || !remove_fixed_params)
                n++;
        }
        for (auto &p : ordinalParameters)
        {
            if (!p.second.fixed())
                p.second.val = p.second.min_val + list[n]*(p.second.max_val - p.second.min_val);
            if (!p.second.fixed() || !remove_fixed_params)
                n++;
        }
        for (auto &p : continuousParameters)
        {
            if (!p.second.fixed())
                p.second.val = p.second.min_val + list[n]*(p.second.max_val - p.second.min_val);
            if (!p.second.fixed() || !remove_fixed_params)
                n++;
        }
    }
    else
    {
        for (auto &p : categorialParameters)
        {
            p.second.val = list[n];
            n++;
        }
        for (auto &p : ordinalParameters)
        {
            p.second.val = list[n];
            n++;
        }
        for (auto &p : continuousParameters)
        {
            p.second.val = list[n];
            n++;
        }
    }
}

float ParameterList::diff(ParameterList &list, bool normalized, bool remove_fixed_params)
{
    if (categorialParameters.size() != list.categorialParameters.size() ||
        ordinalParameters.size() != list.ordinalParameters.size() ||
        continuousParameters.size() != list.continuousParameters.size())
    {
        return 1;
    }

    float dist = 0;
    {
        auto it1 = categorialParameters.begin();
        auto it2 = list.categorialParameters.begin();
        for (int i=0;i<categorialParameters.size();i++)
        {
            if (it1->second.val != it2->second.val)
            {
                dist += 1;
            }
            it1++;
            it2++;
        }
    }
    {
        auto it1 = ordinalParameters.begin();
        auto it2 = list.ordinalParameters.begin();
        for (int i=0;i<ordinalParameters.size();i++)
        {
            if (!it1->second.fixed())
            {
                float sz = it1->second.max_val - it2->second.min_val;
                dist += abs(it2->second.val - it1->second.val)/sz;
            }
            it1++;
            it2++;
        }
    }
    {
        auto it1 = continuousParameters.begin();
        auto it2 = list.continuousParameters.begin();
        for (int i=0;i<continuousParameters.size();i++)
        {
            if (!it1->second.fixed())
            {
                float sz = it1->second.max_val - it2->second.min_val;
                dist += abs(it2->second.val - it1->second.val)/sz;
            }
            it1++;
            it2++;
        }
    }
    return dist / (categorialParameters.size() + ordinalParameters.size() + continuousParameters.size());
}