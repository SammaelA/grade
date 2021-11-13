#include "mygen_parameters.h"


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

        {{ParameterMaskValues::ONE_VALUE,"base_light"}, base_light},
        {{ParameterMaskValues::ONE_VALUE,"base_light_pow"}, base_light_pow},
        
        {{ParameterMaskValues::ONE_VALUE,"dist_power"}, dist_power},
        {{ParameterMaskValues::FULL,"dist_mul"}, dist_mul},
        {{ParameterMaskValues::ONE_VALUE,"base_dist"}, base_dist},
    }; 
    for (auto &p: list)
    {
        if ((double)(p.second.get_max() - p.second.get_min()) < 1e-4 || 
            (double)abs(p.second.get_max() - p.second.get_min()) > 1e10)
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