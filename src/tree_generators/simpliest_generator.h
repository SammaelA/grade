#pragma once
#include "abstract_generator.h"
#include "common_utils/parameter.h"

struct SimpliestTreeStructureParameters : public ParameterSet
{
    int max_depth = 4;
    float4 branch_len = float4(60,20,12,7);
    float4 branch_r = float4(1.75,0.75,0.2, 0.075);
    float4 branch_angle = float4(PI/6, PI/6, PI/4, PI/4);
    float4 branch_count = float4(20,10,8,10);
    float4 branching_start = float4(0.5,0.0,0,0);
    float leaves_count = 0.75;
    float leaf_size = 3;
    virtual ParameterSet *copy() override
    { 
        auto Ps = new SimpliestTreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
    virtual float3 get_tree_max_size() override
    {
        float len = branch_len[0] + branch_len[1] + branch_len[2] + branch_len[3];
        return float3(len,len,len);
    }
    
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
};

class SimpliestTreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void plant_tree(float3 pos, const TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
private:
    void create_tree(Tree *tree, float3 pos, SimpliestTreeStructureParameters &params);
    void create_branch(Tree *tree, Branch *branch, float3 start_pos, float3 base_dir, float3 normal, int level, 
                       SimpliestTreeStructureParameters &params, int &leaves_tries, int &leaves_cnt);
    std::vector<float3> tree_positions;
    std::vector<const TreeTypeData *> types;
};