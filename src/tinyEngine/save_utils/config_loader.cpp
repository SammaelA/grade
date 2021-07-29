#include "parser.h"
#include "ggd_parser.h"
#include "config.h"
#include "../../tree.h"
#include <malloc.h>
extern "C" int parse_file(char *file);

TreeStructureParameters Config::get(std::string name)
{
    auto it = presets.find(name);
    if (it != presets.end())
        return (*it).second;
    else
        return TreeStructureParameters();
}
GroveGenerationData Config::get_ggd(std::string name)
{
    auto it = ggds.find(name);
    if (it != ggds.end())
        return (*it).second;
    else
        return ggds.at("default");
}
bool Config::load_config()
{
    pd.presets_n = 0;
    dat.ggds_c = 0;
    char file[] = "presets.txt";
    parse_file(file);
    for (int i=0;i<pd.presets_n;i++)
    {
        TreeStructureParameters tsp;
        std::string name(pd.presets[i].name);
        Preset &preset = pd.presets[i];
        for (int j=0;j<preset.params_n;j++)
        {
            Param &p = preset.params[j];
            float baseValue;
            float maxValue = -1e10;
            float minValue = 1e10;
            std::vector<float> stateParams;
            int state = -1;
            Normal *normal = nullptr;
            Uniform *uniform = nullptr;
            RandomnessLevel randomnessLevel = NO_RANDOM;
            
            baseValue = p.param_base;
            if (p.param_type > 1)
            {
                for (int k=0;k<p.params;k++)
                {
                    stateParams.push_back(p.param_list[k]);
                }
            }
            if (p.param_type > 2)
            {
                if (p.dist_type == 1)
                {
                    normal = distibutionGenerator.get_normal(p.distr_a,p.distr_b);
                    uniform = distibutionGenerator.get_uniform(0,0);
                }
                else
                {
                    normal = distibutionGenerator.get_normal(0,0);
                    uniform = distibutionGenerator.get_uniform(p.distr_a,p.distr_b);
                    
                }
                switch (p.regen_type)
                {
                case 1:
                    randomnessLevel = NO_RANDOM;
                    break;
                case 2:
                    randomnessLevel = EXPLICIT_REGENERATION;
                    break;
                case 3:
                    randomnessLevel = REGENERATE_ON_STATE_CHANGE;
                    break;
                case 4:
                    randomnessLevel = REGENERATE_ON_GET;
                    break;
                default:
                    break;
                }
            }
            if (p.param_type > 3)
            {
                minValue = p.param_from;
                maxValue = p.param_to;
            }
            Parameter<float> res = Parameter<float>(baseValue,stateParams,randomnessLevel,normal,minValue,maxValue);
            if (p.dist_type != 1)
                res = Parameter<float>(baseValue,stateParams,randomnessLevel,uniform,minValue,maxValue);
            std::string nm = p.name;
            if (nm == "max_depth")
            {
                tsp.max_depth.set_no_override_minmax(res);
            }
            else if (nm == "max_segments")
            {
                tsp.max_segments.set_no_override_minmax(res);
            }
            else if (nm == "max_branching")
            {
                tsp.max_branching.set_no_override_minmax(res);
            }
            else if (nm == "growth_iterations")
            {
                tsp.growth_iterations.set_no_override_minmax(res);;
            }

            else if (nm == "scale")
            {
                tsp.scale.set_no_override_minmax(res);
            }
            else if (nm == "seg_len_mult")
            {
                tsp.seg_len_mult.set_no_override_minmax(res);
            }
            else if (nm == "leaf_size_mult")
            {
                tsp.leaf_size_mult.set_no_override_minmax(res);
            }
            else if (nm == "base_r")
            {
                tsp.base_r.set_no_override_minmax(res);
            }
            else if (nm == "r_split_save_pow")
            {
                tsp.r_split_save_pow.set_no_override_minmax(res);
            }

            else if (nm == "dir_conserv")
            {
                tsp.dir_conserv.set_no_override_minmax(res);
            }
            else if (nm == "plane_conserv")
            {
                tsp.plane_conserv.set_no_override_minmax(res);
            }
            else if (nm == "spread")
            {
                tsp.spread.set_no_override_minmax(res);
            }
            else if (nm == "phototrop")
            {
                tsp.phototrop.set_no_override_minmax(res);
            }
            else if (nm == "gravitrop")
            {
                tsp.gravitrop.set_no_override_minmax(res);
            }
            else if (nm == "dir_random")
            {
                tsp.dir_random.set_no_override_minmax(res);
            }
            else if (nm == "base_angle")
            {
                tsp.base_angle.set_no_override_minmax(res);
            }
            else if (nm == "base_angle_q")
            {
                tsp.base_angle_q.set_no_override_minmax(res);
            }

            else if (nm == "seg_dir_conserv")
            {
                tsp.seg_dir_conserv.set_no_override_minmax(res);
            }
            else if (nm == "seg_plane_conserv")
            {
                tsp.seg_plane_conserv.set_no_override_minmax(res);
            }
            else if (nm == "seg_spread")
            {
                tsp.seg_spread.set_no_override_minmax(res);
            }
            else if (nm == "seg_phototrop")
            {
                tsp.seg_phototrop.set_no_override_minmax(res);
            }
            else if (nm == "seg_gravitrop")
            {
                tsp.seg_gravitrop.set_no_override_minmax(res);
            }
            else if (nm == "seg_dir_random")
            {
                tsp.seg_dir_random.set_no_override_minmax(res);
            }
            else if (nm == "seg_bend")
            {
                tsp.seg_bend.set_no_override_minmax(res);
            }
            else if (nm == "seg_bend_pow")
            {
                tsp.seg_bend_pow.set_no_override_minmax(res);
            }

            else if (nm == "base_branch_feed")
            {
                tsp.base_branch_feed.set_no_override_minmax(res);
            }
            else if (nm == "base_seg_feed")
            {
                tsp.base_seg_feed.set_no_override_minmax(res);
            }
            else if (nm == "feed_distribution_min_weight")
            {
                tsp.feed_distribution_min_weight.set_no_override_minmax(res);
            }
            else if (nm == "feed_distribution_d_weight")
            {
                tsp.feed_distribution_d_weight.set_no_override_minmax(res);
            }
            else if (nm == "top_growth_bonus")
            {
                tsp.top_growth_bonus.set_no_override_minmax(res);
            }

            else if (nm == "light_precision")
            {
                tsp.light_precision.set_no_override_minmax(res);
            }
            else if (nm == "branch_removal")
            {
                tsp.branch_removal.set_no_override_minmax(res);
            }
            else if (nm == "branch_grow_decrease_q")
            {
                tsp.branch_grow_decrease_q.set_no_override_minmax(res);
            }
            else if (nm == "segment_grow_decrease_q")
            {
                tsp.segment_grow_decrease_q.set_no_override_minmax(res);
            }

            else if (nm == "min_branching_chance")
            {
                tsp.min_branching_chance.set_no_override_minmax(res);
            }
            else if (nm == "max_branching_chance")
            {
                tsp.max_branching_chance.set_no_override_minmax(res);
            }
            else if (nm == "branching_power")
            {
                tsp.branching_power.set_no_override_minmax(res);
            }

            else if (nm == "r_deformation_levels")
            {
                tsp.r_deformation_levels.set_no_override_minmax(res);
            }
            else if (nm == "r_deformation_points")
            {
                tsp.r_deformation_points.set_no_override_minmax(res);
            }
            else if (nm == "r_deformation_power")
            {
                tsp.r_deformation_power.set_no_override_minmax(res);
            }
            else if (nm == "max_segments")
            {
                tsp.max_segments.set_no_override_minmax(res);
            }
            else if (nm == "seg_len_mult")
            {
                tsp.seg_len_mult.set_no_override_minmax(res);
            }
            else if (nm == "base_seg_feed")
            {
                tsp.base_seg_feed.set_no_override_minmax(res);
            }
            else if (nm == "base_branch_feed")
            {
                tsp.base_branch_feed.set_no_override_minmax(res);
            }
            else if (nm == "min_branching_chance")
            {
                tsp.min_branching_chance.set_no_override_minmax(res);
            }
            free(p.name);
            free(p.param_list);
        }
        presets.emplace(name,tsp);
    }
}