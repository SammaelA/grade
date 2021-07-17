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
            Distribution *randomizer = nullptr;
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
                    randomizer = distibutionGenerator.get_normal(p.distr_a,p.distr_b);
                else
                    randomizer = distibutionGenerator.get_uniform(p.distr_a,p.distr_b);
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
            Parameter<float> res = Parameter<float>(baseValue,stateParams,randomnessLevel,randomizer,minValue,maxValue);
            std::string nm = p.name;
            if (nm == "max_depth")
            {
                tsp.max_depth = TreeStructureParameters::from_float(res);
            }
            else if (nm == "max_segments")
            {
                tsp.max_segments = TreeStructureParameters::from_float(res);
            }
            else if (nm == "max_branching")
            {
                tsp.max_branching = TreeStructureParameters::from_float(res);
            }
            else if (nm == "growth_iterations")
            {
                tsp.growth_iterations = TreeStructureParameters::from_float(res);
            }

            else if (nm == "scale")
            {
                tsp.scale = res;
            }
            else if (nm == "seg_len_mult")
            {
                tsp.seg_len_mult = res;
            }
            else if (nm == "leaf_size_mult")
            {
                tsp.leaf_size_mult = res;
            }
            else if (nm == "base_r")
            {
                tsp.base_r = res;
            }
            else if (nm == "r_split_save_pow")
            {
                tsp.r_split_save_pow = res;
            }

            else if (nm == "dir_conserv")
            {
                tsp.dir_conserv = res;
            }
            else if (nm == "plane_conserv")
            {
                tsp.plane_conserv = res;
            }
            else if (nm == "spread")
            {
                tsp.spread = res;
            }
            else if (nm == "phototrop")
            {
                tsp.phototrop = res;
            }
            else if (nm == "gravitrop")
            {
                tsp.gravitrop = res;
            }
            else if (nm == "dir_random")
            {
                tsp.dir_random = res;
            }
            else if (nm == "base_angle")
            {
                tsp.base_angle = res;
            }
            else if (nm == "base_angle_q")
            {
                tsp.base_angle_q = res;
            }

            else if (nm == "seg_dir_conserv")
            {
                tsp.seg_dir_conserv = res;
            }
            else if (nm == "seg_plane_conserv")
            {
                tsp.seg_plane_conserv = res;
            }
            else if (nm == "seg_spread")
            {
                tsp.seg_spread = res;
            }
            else if (nm == "seg_phototrop")
            {
                tsp.seg_phototrop = res;
            }
            else if (nm == "seg_gravitrop")
            {
                tsp.seg_gravitrop = res;
            }
            else if (nm == "seg_dir_random")
            {
                tsp.seg_dir_random = res;
            }
            else if (nm == "seg_bend")
            {
                tsp.seg_bend = res;
            }
            else if (nm == "seg_bend_pow")
            {
                tsp.seg_bend_pow = res;
            }

            else if (nm == "base_branch_feed")
            {
                tsp.base_branch_feed = res;
            }
            else if (nm == "base_seg_feed")
            {
                tsp.base_seg_feed = res;
            }
            else if (nm == "feed_distribution_min_weight")
            {
                tsp.feed_distribution_min_weight = res;
            }
            else if (nm == "feed_distribution_d_weight")
            {
                tsp.feed_distribution_d_weight = res;
            }
            else if (nm == "top_growth_bonus")
            {
                tsp.top_growth_bonus = res;
            }

            else if (nm == "light_precision")
            {
                tsp.light_precision = res;
            }
            else if (nm == "branch_removal")
            {
                tsp.branch_removal = res;
            }
            else if (nm == "branch_grow_decrease_q")
            {
                tsp.branch_grow_decrease_q = res;
            }
            else if (nm == "segment_grow_decrease_q")
            {
                tsp.segment_grow_decrease_q = res;
            }

            else if (nm == "min_branching_chance")
            {
                tsp.min_branching_chance = res;
            }
            else if (nm == "max_branching_chance")
            {
                tsp.max_branching_chance = res;
            }
            else if (nm == "branching_power")
            {
                tsp.branching_power = res;
            }

            else if (nm == "r_deformation_levels")
            {
                tsp.r_deformation_levels = TreeStructureParameters::from_float(res);
            }
            else if (nm == "r_deformation_points")
            {
                tsp.r_deformation_points = TreeStructureParameters::from_float(res);
            }
            else if (nm == "r_deformation_power")
            {
                tsp.r_deformation_power = res;
            }
            else if (nm == "max_segments")
            {
                tsp.max_segments = TreeStructureParameters::from_float(res);
            }
            else if (nm == "seg_len_mult")
            {
                tsp.seg_len_mult = res;
            }
            else if (nm == "base_seg_feed")
            {
                tsp.base_seg_feed = res;
            }
            else if (nm == "base_branch_feed")
            {
                tsp.base_branch_feed = res;
            }
            else if (nm == "min_branching_chance")
            {
                tsp.min_branching_chance = res;
            }
            free(p.name);
            free(p.param_list);
        }
        presets.emplace(name,tsp);
    }
}