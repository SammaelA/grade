#pragma once

#include "graphics_utils/volumetric_occlusion.h"
#include "abstract_generator.h"
#include "common_utils/parameter.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "common_utils/octree.h"
#include <vector>
#include <list>
#include <atomic>

struct GETreeParameters : public ParametersSet
{
    float lambda = 0.52;
    float k = 0.75;//part of joint that can create child branches
    int tau = 6;
    float ro = 1.0;
    float X0 = 2;
    float Xm = 100;
    float r = 0.37;
    int alpha = 4;
    float sigma = 0.5;
    float mu = 0.75;
    float b_min = 1.8;
    float b_max = 2.2;
    float r_s = 0.05;
    float rs_size_factor = 0.02;
    int remove_min_level = 1;

    float base_r = 0.035;
    int max_branches = 1;
    int occlusion_pyramid_d = 7;
    float r_pow = 2.2;
    int sp_points_base = 4;
    float branching_angle_min = 0;
    float branching_angle_max = PI/3;
    int max_iterations = 100;
    float leaf_size_mult = 3.5;
    float leaves_cnt = 1.0;
    int max_joints_in_branch = 32;
    float resource_mult = 6.0;
    float res_q = 0.5;
    float leaves_max_r = 5;//if radius in node > leaves_max_r*base_r leaf will not be created on this node
    float leaves_angle_a = 0.3;
    float leaves_angle_b = 0.4;
    int root_type = 1;
    glm::vec4 tropism_params = glm::vec4(1,1,1.0/15,1.5);
    glm::vec2 tropism_min_max = glm::vec2(-5,5);
    float branching_tropims_mult = 0.67;
    int tropism_level_base = 1;
    float res_decrease_step = 0.09;
    float res_decrease_min = 0.03;

    virtual glm::vec3 get_tree_max_size() override
    {
        if (root_type == 0 || root_type == 2)
            return ro*glm::vec3(1.5*Xm, Xm + 30, 1.5*Xm);
        else if (root_type == 1)
            return ro*glm::vec3(0.6*Xm, 1.25*Xm + 30, 0.6*Xm);
        else 
            return ro*glm::vec3(2.0f*Xm, 1.5f*Xm, 2.0f*Xm);
    }
    virtual float get_scale_factor() override 
    {
        return ro;
    }
    virtual ParametersSet *copy() override
    { 
        auto Ps = new GETreeParameters();
        *Ps = *this;
        return Ps;
    };
    virtual void save_to_blk(Block &b) override;
    virtual void load_from_blk(Block &b) override;
    virtual void RW_parameter_list(bool write, ParameterList &list) override;
};