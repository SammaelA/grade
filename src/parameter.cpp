#include "parameter.h"
#include "tinyEngine/utility.h"


Parameter<int> TreeStructureParameters::from_float(Parameter<float> source)
{
    std::vector<float> qs = source.state_qs;
    
    Parameter<int> par((int)source.baseValue,(int)source.minValue,(int)source.maxValue, qs, source.a, source.sigma, 
                       source.from, source.to, source.normal_part, source.randomnessLevel,
                       dynamic_cast<Normal*>(source.normal), dynamic_cast<Uniform*>(source.uniform));
    return par;
}
Parameter<float> TreeStructureParameters::from_int(Parameter<int> source)
{
    std::vector<float> qs = source.state_qs;
    
    Parameter<float> par((float)source.baseValue,(float)source.minValue,(float)source.maxValue, qs, source.a, source.sigma, 
                       source.from, source.to, source.normal_part, source.randomnessLevel,
                       dynamic_cast<Normal*>(source.normal), dynamic_cast<Uniform*>(source.uniform));
    return par;
}

void TreeStructureParameters::get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                                                 ParameterVariablesSet v_set)
{
    list = {
        {{ParameterMaskValues::CONSTANT,"max_depth"}, max_depth},
        {{ParameterMaskValues::CONSTANT,"max_segments"}, max_segments},
        {{ParameterMaskValues::CONSTANT,"max_branching"}, max_branching},
        {{ParameterMaskValues::CONSTANT,"growth_iterations"}, growth_iterations},

        {{ParameterMaskValues::CONSTANT,"scale"}, scale},
        {{ParameterMaskValues::FULL,"seg_len_mult"}, seg_len_mult},
        {{ParameterMaskValues::ONE_VALUE,"leaf_size_mult"}, leaf_size_mult},
        {{ParameterMaskValues::CONSTANT,"base_r"}, base_r},
        {{ParameterMaskValues::LIST_OF_VALUES,"r_split_save_pow"}, r_split_save_pow},

        {{ParameterMaskValues::FULL,"dir_conserv"}, dir_conserv},
        {{ParameterMaskValues::FULL,"plane_conserv"}, plane_conserv},
        {{ParameterMaskValues::FULL,"spread"}, spread},
        {{ParameterMaskValues::FULL,"phototrop"}, phototrop},
        {{ParameterMaskValues::FULL,"gravitrop"}, gravitrop},
        {{ParameterMaskValues::FULL,"dir_random"}, dir_random},
        {{ParameterMaskValues::ONE_VALUE,"base_angle"}, base_angle},
        {{ParameterMaskValues::FULL,"base_angle_q"}, base_angle_q},

        {{ParameterMaskValues::FULL,"seg_dir_conserv"}, seg_dir_conserv},
        {{ParameterMaskValues::FULL,"seg_plane_conserv"}, seg_plane_conserv},
        {{ParameterMaskValues::FULL,"seg_spread"}, seg_spread},
        {{ParameterMaskValues::FULL,"seg_phototrop"}, seg_phototrop},
        {{ParameterMaskValues::FULL,"seg_gravitrop"}, seg_gravitrop},
        {{ParameterMaskValues::FULL,"seg_dir_random"}, seg_dir_random},
        {{ParameterMaskValues::FULL,"seg_bend"}, seg_bend},
        {{ParameterMaskValues::FULL,"seg_bend_pow"}, seg_bend_pow},
        
        {{ParameterMaskValues::FULL,"base_branch_feed"}, base_branch_feed},
        {{ParameterMaskValues::FULL,"base_seg_feed"}, base_seg_feed},
        {{ParameterMaskValues::ONE_VALUE,"feed_distribution_min_weight"}, feed_distribution_min_weight},
        {{ParameterMaskValues::ONE_VALUE,"feed_distribution_d_weight"}, feed_distribution_d_weight},
        {{ParameterMaskValues::ONE_VALUE,"top_growth_bonus"}, top_growth_bonus},

        {{ParameterMaskValues::ONE_VALUE,"light_precision"}, light_precision},
        {{ParameterMaskValues::FULL,"branch_removal"}, branch_removal},
        {{ParameterMaskValues::ONE_VALUE,"branch_grow_decrease_q"}, branch_grow_decrease_q},
        {{ParameterMaskValues::ONE_VALUE,"segment_grow_decrease_q"}, segment_grow_decrease_q},

        {{ParameterMaskValues::LIST_OF_VALUES,"min_branching_chance"}, min_branching_chance},
        {{ParameterMaskValues::ONE_VALUE,"max_branching_chance"}, max_branching_chance},
        {{ParameterMaskValues::LIST_OF_VALUES,"branching_power"}, branching_power},

        {{ParameterMaskValues::ONE_VALUE,"r_deformation_levels"}, r_deformation_levels},
        {{ParameterMaskValues::CONSTANT,"r_deformation_points"}, r_deformation_points},
        {{ParameterMaskValues::CONSTANT,"r_deformation_power"}, r_deformation_power},
    }; 
    for (auto &p: list)
    {
        if ((double)(p.second.maxValue - p.second.minValue) < 1e-4 || 
            (double)abs(p.second.maxValue - p.second.minValue) > 1e10)
        {
            p.first.val = ParameterMaskValues::CONSTANT;
        }
        if (v_set == ParameterVariablesSet::ONLY_BASE_VALUES && p.first.val != ParameterMaskValues::CONSTANT)
        {
            p.first.val = ParameterMaskValues::ONE_VALUE;
        }
        else if (v_set == ParameterVariablesSet::BASE_VALUES_AND_QS && p.first.val == ParameterMaskValues::FULL)
        {
            p.first.val = ParameterMaskValues::LIST_OF_VALUES;
        }
        /*
        if (p.first.val == ParameterMaskValues::CONSTANT)
        {
            logerr("CONSTANT");
        }
        if (p.first.val == ParameterMaskValues::FULL)
        {
            logerr("FULL");
        }
        if (p.first.val == ParameterMaskValues::LIST_OF_VALUES)
        {
            logerr("LIST_OF_VALUES");
        }
        if (p.first.val == ParameterMaskValues::ONE_VALUE)
        {
            logerr("ONE_VALUE");
        }
        logerr("%s %s %d \n",p.first.name.c_str(), p.second.to_string().c_str(),
               (int)p.second.randomnessLevel);*/
    } 
}
void TreeStructureParameters::get_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, ParameterVariablesSet v_set)
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
void TreeStructureParameters::load_from_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data, ParameterVariablesSet v_set)
{
    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    get_parameter_list(list,v_set);
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