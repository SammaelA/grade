#pragma once
#include "common_utils/parameter.h"

struct TreeStructureParameters : public ParametersSet
{
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

    Parameter<float> base_light;
    Parameter<float> base_light_pow;

    Parameter<float> dist_power;
    Parameter<float> dist_mul;
    Parameter<float> base_dist;
    virtual void set_state(int state) override
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

        base_light.set_state(state);
        base_light_pow.set_state(state);
        
        dist_power.set_state(state);
        dist_mul.set_state(state);
        base_dist.set_state(state);
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
                                r_deformation_power(0, std::vector<float>{0,0}, REGENERATE_ON_GET, distibutionGenerator.get_normal(0,0.025), 0, 1),
                        
                                base_light(40,10,150),
                                base_light_pow(1,0.5,2.5),
                            
                                dist_power(1,0.5,2.5),
                                dist_mul(0.4,std::vector<float>{0, 0.3, 0.5, 0.4},REGENERATE_ON_GET, distibutionGenerator.get_normal(0,0.25),0,4),
                                base_dist(40,0,150)
    {
    }
    virtual void get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                            ParameterVariablesSet v_set = ParameterVariablesSet::ALL_VALUES) override;
    virtual glm::vec3 get_tree_max_size() override
    {
        set_state(0);
        return 2.0f*glm::vec3(seg_len_mult()*max_segments());
    }
    virtual ParametersSet *copy() override
    {
        auto Ps = new TreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
};
template<typename T>
class BaseParameterSetWrapper
{
    std::vector<T> params;
    int state = 0;
    public:
    BaseParameterSetWrapper()
    {
        params.push_back(T());
    }
    BaseParameterSetWrapper(T base, int n_params)
    {
        for (int i=0;i<n_params;i++)
        {
            params.push_back(base);
            params.back().set_state(i);
        }
    }
    void set_state(int _state) {state = _state;}
    T &operator()() {return params[state];}
};

using ParameterSetWrapper = BaseParameterSetWrapper<TreeStructureParameters>;