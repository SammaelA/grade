#include "parameter.h"
#include "tinyEngine/utility.h"


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

void ParametersSet::get_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, ParameterVariablesSet v_set)
{
    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    get_parameter_list(list,v_set);
    for (auto &p: list)
    {
        ParameterDesc desc(p.second);
        desc.mask = p.first.val;
        desc.name = p.first.name;
        desc.randomnessLevel = p.second.randomnessLevel;
        desc.var_count = 0;

        if (p.first.val != ParameterMaskValues::CONSTANT)
        {
            desc.var_count += 1;
            float val = p.second.baseValue;
            if ((double)(p.second.maxValue - p.second.minValue) < 1e-4)
            {
                logerr("max val should be bigger than min val for non-constant parameters. Min val will be ignored %f %f",
                p.second.maxValue, p.second.minValue);
                val = 1;
            }
            else
            {
                val = (val - p.second.minValue) / (p.second.maxValue - p.second.minValue);
            }
            data.push_back(val);
            
            if (p.first.val == ParameterMaskValues::LIST_OF_VALUES)
            {
                desc.var_count += p.second.state_qs.size();
                for (float &q : p.second.state_qs)
                    data.push_back(q);
            }
            else if (p.first.val == ParameterMaskValues::FULL)
            {
                desc.var_count += p.second.state_qs.size();
                for (float &q : p.second.state_qs)
                {
                    data.push_back(q);
                    //if (desc.name == "base_seg_feed")
                    //logerr("qs = %f",q);
                }
                desc.var_count += 4;
                
                //data.push_back(p.second.a);
                float delta = (p.second.maxValue - p.second.minValue);
                data.push_back(p.second.sigma/delta);
                float from_n = CLAMP((p.second.from + delta) / (2*delta), 0, 1);
                float delta_n = CLAMP((p.second.to + delta) / (2*delta), 0, 1);
                data.push_back(from_n);
                data.push_back(delta_n);
                data.push_back(p.second.normal_part);
            }
        }

        mask.push_back(desc);
    }
}
void ParametersSet::load_from_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, ParameterVariablesSet v_set)
{
    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    this->get_parameter_list(list,v_set);
    if (mask.size() != list.size())
    {
        logerr("unable to load tree structure parameters. Mask size %d list size %d",mask.size(), list.size());
        return;
    }
    int cur_pos = 0;
    for (int i=0;i<mask.size();i++)
    {
        if (cur_pos + mask[i].var_count > data.size())
        {
            logerr("unable to load tree structure parameters. Data size %d, should be %d",
                   data.size(), cur_pos + mask[i].var_count);
            return;
        }
        float minValue = mask[i].original.minValue;
        float maxValue = mask[i].original.maxValue;
        float val = mask[i].original.baseValue;
        std::vector<float> state_qs = mask[i].original.state_qs;
        float a = mask[i].original.a;
        float sigma = mask[i].original.sigma;
        float from = mask[i].original.from;
        float to = mask[i].original.to;
        float normal_part = mask[i].original.normal_part;
        RandomnessLevel rl = mask[i].randomnessLevel;
        if (mask[i].mask == ParameterMaskValues::CONSTANT)
            continue;
        else if (mask[i].mask == ParameterMaskValues::ONE_VALUE)
        {
            val = data[cur_pos]*(maxValue - minValue) + minValue;
            list[i].second = Parameter<float>(val,minValue,maxValue);
        }
        else if (mask[i].mask == ParameterMaskValues::LIST_OF_VALUES)
        {
            val = data[cur_pos]*(maxValue - minValue) + minValue;
            std::vector<float> state_vals;
            state_qs = {};
            for (int j=1; j<mask[i].var_count; j++)
            {
                state_qs.push_back(data[cur_pos+j]);
            }
        }
        else if (mask[i].mask == ParameterMaskValues::FULL)
        {
            val = data[cur_pos]*(maxValue - minValue) + minValue;
            state_qs = {};
            for (int j=1; j<mask[i].var_count-4; j++)
            {
                state_qs.push_back(data[cur_pos+j]);
                    //if (mask[i].name == "base_seg_feed")
                    //logerr("qs = %f %d",data[cur_pos+j], cur_pos+j);
            }

            int p = cur_pos + mask[i].var_count-4;

            a = 0;//data[p];
            float delta = (maxValue - minValue);
            sigma = data[p]*delta;
            float from_n = data[p+1];
            float delta_n = data[p+2];
            normal_part = data[p+3];  

            from = 2*delta*from_n - delta;
            to = 2*delta*delta_n - delta;
        }
        list[i].second = Parameter<float>(val,minValue,maxValue,state_qs,a,sigma,from,to,
                                          normal_part,rl,nullptr,nullptr);
        cur_pos += mask[i].var_count;        
    }
    //get_parameter_list(list,v_set);
}