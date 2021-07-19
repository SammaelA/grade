#include "parameter.h"
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
void TreeStructureParameters::get_mask_and_data(std::vector<ParameterMaskValues> &mask, std::vector<double> &data)
{

}