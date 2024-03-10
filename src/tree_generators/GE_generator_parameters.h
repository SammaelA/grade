#pragma once

#include "graphics_utils/volumetric_occlusion.h"
#include "abstract_generator.h"
#include "common_utils/parameter.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "common_utils/octree.h"
#include <vector>
#include <list>
#include <atomic>

struct GETreeParameters : public ParameterSet
{
    float lambda = 0.52;
    float k = 0.75;//part of joint that can create child branches
    int tau = 6;//unused
    float ro = 1.0;
    float X0 = 2;
    float Xm = 100;
    float r = 0.37;
    int alpha = 4;//unused
    float sigma = 0.5;//unused
    float mu = 0.75;
    float b_min = 1.8;//unused
    float b_max = 2.2;//unused
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
    float4 tropism_params = float4(1,1,1.0/15,1.5);
    float2 tropism_min_max = float2(-5,5);
    float branching_tropims_mult = 0.67;
    int tropism_level_base = 0;
    float res_decrease_step = 0.09;
    float res_decrease_min = 0.03;
    float2 initial_trunk_scale = float2(1,1);//horizontal, vertical scale
    float trunk_bonus_radius = 0.1;
    float trunk_bonus_radius_mod = 1;
    float2 distance_fine_start = float2(1000,1000);
    float2 distance_fine_slope = float2(1000,1000);

    //Beta test
    float L0 = 30;
    float R0 = 4;
    int B_cnt = 8;
    float dX = 1;
    float4 iR = float4(0.67, 0.67, 0.67, 0.67);
    float4 dR = float4(0.85, 0.85, 0.85, 0.85);
    float4 Lt = float4(1,3,5,7);
    float4 Lb = float4(5,10,15,20);
    float4 phi = float4(0.0,0.5,1,1.5);// phi/(2*PI)
    float4 psi = float4(0.1,0.3,0.5,0.7);//psi /(PI/2)

    virtual float3 get_tree_max_size() override
    {
        if (root_type == 0 || root_type == 2)
            return ro*float3(1.5*Xm, Xm + 30, 1.5*Xm);
        else if (root_type == 1)
            return ro*float3(0.6*Xm, 1.25*Xm + 0.4*Xm*initial_trunk_scale.y, 0.6*Xm);
        else if (root_type == 3)
        {
            float max_y = MAX(0, L0) + B_cnt*MAX(MAX(Lt.x, Lt.y), MAX(Lt.z, Lt.w));
            return ro*float3(0.6*Xm, 1.25*Xm + max_y*initial_trunk_scale.y, 0.6*Xm);
        }
        else 
            return ro*float3(2.0f*Xm, 1.5f*Xm, 2.0f*Xm);
    }
    virtual float get_scale_factor() override 
    {
        return ro;
    }
    virtual ParameterSet *copy() override
    { 
        auto Ps = new GETreeParameters();
        *Ps = *this;
        return Ps;
    };
    
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
};