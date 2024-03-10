#pragma once
#include "abstract_generator.h"
#include "common_utils/parameter.h"
struct SimpleTreeStructureParameters : public ParameterSet
{
    int max_depth = 4;
    float4 segment_size = float4(3.0,1.5,0.5,0.5);
    float4 segment_count = float4(15, 12, 10, 15);
    float4 base_thickness = float4(2,1.4,0.6,0.4);
    float base_dir_mult = 5;
    float rand_dir_mult = 1;
    float    up_dir_mult = 1;
    float    down_dir_mult = 1;
    float    up_dir_pow = 1;
    float    down_dir_pow = 2;
    float    sectors_count = 6;
    float    base_branch_angle = PI/4;
    float    branching_chance = 0.75;
    float    leaves_chance = 1;
    float    leaf_mult = 2;

    
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
    virtual float3 get_tree_max_size() override
    {
      return 2.0f*float3(segment_size[0]*segment_count[0]);
    }
    virtual ParameterSet *copy() override
    { 
        auto Ps = new SimpleTreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
};

class SimpleTreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void plant_tree(float3 pos, const TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
private:
    void create_tree(Tree *tree, float3 pos);
    void create_branch(Tree *tree, Branch *branch, float3 start_pos, float3 base_dir, float3 normal, int level, 
                       float base_r, float leaves_chance);
    SimpleTreeStructureParameters params;
    std::vector<float3> tree_positions;
    std::vector<const TreeTypeData *> types;
};