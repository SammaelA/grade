#pragma once
#include "abstract_generator.h"
#include "common_utils/parameter.h"

struct SimpliestTreeStructureParameters : public ParameterSet
{
    int max_depth = 4;
    glm::vec4 branch_len = glm::vec4(60,20,12,7);
    glm::vec4 branch_r = glm::vec4(1.75,0.75,0.2, 0.075);
    glm::vec4 branch_angle = glm::vec4(PI/6, PI/6, PI/4, PI/4);
    glm::vec4 branch_count = glm::vec4(20,10,8,10);
    glm::vec4 branching_start = glm::vec4(0.5,0.0,0,0);
    float leaves_count = 0.75;
    float leaf_size = 3;
    virtual ParameterSet *copy() override
    { 
        auto Ps = new SimpliestTreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
    virtual glm::vec3 get_tree_max_size() override
    {
        float len = branch_len[0] + branch_len[1] + branch_len[2] + branch_len[3];
        return glm::vec3(len,len,len);
    }
    
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
};

class SimpliestTreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void plant_tree(glm::vec3 pos, const TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
private:
    void create_tree(Tree *tree, glm::vec3 pos, SimpliestTreeStructureParameters &params);
    void create_branch(Tree *tree, Branch *branch, glm::vec3 start_pos, glm::vec3 base_dir, glm::vec3 normal, int level, 
                       SimpliestTreeStructureParameters &params, int &leaves_tries, int &leaves_cnt);
    std::vector<glm::vec3> tree_positions;
    std::vector<const TreeTypeData *> types;
};