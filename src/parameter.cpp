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
void TreeStructureParameters::get_parameter_list(std::vector<std::pair<ParameterMaskValues,Parameter<float> &>> &list)
{
    list = {
        {ParameterMaskValues::CONSTANT, max_depth},
        {ParameterMaskValues::FULL, max_segments},
        {ParameterMaskValues::CONSTANT, max_branching},
        {ParameterMaskValues::CONSTANT, growth_iterations},

        {ParameterMaskValues::CONSTANT, scale},
        {ParameterMaskValues::FULL, seg_len_mult},
        {ParameterMaskValues::ONE_VALUE, leaf_size_mult},
        {ParameterMaskValues::FULL, base_r},
        {ParameterMaskValues::LIST_OF_VALUES, r_split_save_pow},

        {ParameterMaskValues::FULL, dir_conserv},
        {ParameterMaskValues::FULL, plane_conserv},
        {ParameterMaskValues::FULL, spread},
        {ParameterMaskValues::FULL, phototrop},
        {ParameterMaskValues::FULL, gravitrop},
        {ParameterMaskValues::FULL, dir_random},
        {ParameterMaskValues::ONE_VALUE, base_angle},
        {ParameterMaskValues::FULL, base_angle_q},

        {ParameterMaskValues::FULL, seg_dir_conserv},
        {ParameterMaskValues::FULL, seg_plane_conserv},
        {ParameterMaskValues::FULL, seg_spread},
        {ParameterMaskValues::FULL, seg_phototrop},
        {ParameterMaskValues::FULL, seg_gravitrop},
        {ParameterMaskValues::FULL, seg_dir_random},
        {ParameterMaskValues::FULL, seg_bend},
        {ParameterMaskValues::FULL, seg_bend_pow},
        
        {ParameterMaskValues::FULL, base_branch_feed},
        {ParameterMaskValues::FULL, base_seg_feed},
        {ParameterMaskValues::ONE_VALUE, feed_distribution_min_weight},
        {ParameterMaskValues::ONE_VALUE, feed_distribution_d_weight},
        {ParameterMaskValues::ONE_VALUE, top_growth_bonus},

        {ParameterMaskValues::ONE_VALUE, light_precision},
        {ParameterMaskValues::FULL, branch_removal},
        {ParameterMaskValues::ONE_VALUE, branch_grow_decrease_q},
        {ParameterMaskValues::ONE_VALUE, segment_grow_decrease_q},

        {ParameterMaskValues::LIST_OF_VALUES, min_branching_chance},
        {ParameterMaskValues::ONE_VALUE, max_branching_chance},
        {ParameterMaskValues::LIST_OF_VALUES, branching_power},

        {ParameterMaskValues::ONE_VALUE, r_deformation_levels},
        {ParameterMaskValues::LIST_OF_VALUES, r_deformation_points},
        {ParameterMaskValues::LIST_OF_VALUES, r_deformation_power},
    };  
}
void TreeStructureParameters::get_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data)
{
    std::vector<std::pair<ParameterMaskValues,Parameter<float> &>> list;
    get_parameter_list(list);
    for (auto &p: list)
    {
        ParameterDesc desc;
        desc.mask = p.first;
        desc.minValue = p.second.minValue;
        desc.maxValue = p.second.maxValue;
        desc.count = 0;
        if (p.first != ParameterMaskValues::CONSTANT)
        {
            desc.count += 1;
            data.push_back(p.second.baseValue);
            
            if (p.first == ParameterMaskValues::LIST_OF_VALUES)
            {
                desc.count += p.second.state_qs.size();
                for (float &q : p.second.state_qs)
                    data.push_back(q);
            }
            else if (p.first == ParameterMaskValues::FULL)
            {
                desc.count += p.second.state_qs.size();
                for (float &q : p.second.state_qs)
                    data.push_back(q);

                desc.count += 5;
                
                data.push_back(p.second.a);
                data.push_back(p.second.sigma);
                data.push_back(p.second.from);
                data.push_back(p.second.to);
                data.push_back(p.second.normal_part);
            }
        }

        mask.push_back(desc);
    }
}
void TreeStructureParameters::load_from_mask_and_data(std::vector<ParameterDesc> &mask, std::vector<double> &data)
{
    std::vector<std::pair<ParameterMaskValues,Parameter<float> &>> list;
    get_parameter_list(list);
    if (mask.size() != list.size())
    {
        logerr("unable to load tree structure parameters. Mask size %d list size %d",mask.size(), list.size());
        return;
    }
    int cur_pos = 0;
    for (int i=0;i<mask.size();i++)
    {
        if (cur_pos + mask[i].count > data.size())
        {
            logerr("unable to load tree structure parameters. Data size %d, should be %d",
                   data.size(), cur_pos + mask[i].count);
            return;
        }
        if (mask[i].mask == ParameterMaskValues::CONSTANT)
            continue;
        else if (mask[i].mask == ParameterMaskValues::ONE_VALUE)
        {
            float val = data[cur_pos];
            list[i].second = Parameter<float>(val,mask[i].minValue,mask[i].maxValue);
        }
        else if (mask[i].mask == ParameterMaskValues::LIST_OF_VALUES)
        {
            float val = data[cur_pos];
            std::vector<float> state_vals;
            for (int j=1; j<mask[i].count; j++)
            {
                state_vals.push_back(val*data[cur_pos+j]);
            }
            list[i].second = Parameter<float>(val,state_vals, mask[i].minValue,mask[i].maxValue);
        }
        else if (mask[i].mask == ParameterMaskValues::FULL)
        {
            float val = data[cur_pos];
            std::vector<float> state_qs;
            for (int j=1; j<mask[i].count-5; j++)
            {
                state_qs.push_back(data[cur_pos+j]);
            }
            int p = cur_pos + mask[i].count-5;

            float a = data[p];
            float sigma = data[p+1];
            float from = data[p+2];
            float to = data[p+3];
            float normal_part = data[p+4];

            list[i].second = Parameter<float>(val,mask[i].minValue,mask[i].maxValue,state_qs,a,sigma,from,to,
                                              normal_part,RandomnessLevel::REGENERATE_ON_GET,nullptr,nullptr);
            //if (state_qs.size() >=3 )
            //logerr("%f (%f %f %f) %f %f %f %f %f",val,state_qs[0], state_qs[1], state_qs[2],a,sigma,from,to,normal_part);
            //else 
            //logerr("%f (%d) %f %f %f %f %f",val,state_qs.size(),a,sigma,from,to,normal_part);
        }
        cur_pos += mask[i].count;
    }
}