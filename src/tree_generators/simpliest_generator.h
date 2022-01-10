#pragma once
#include "abstract_generator.h"
#include "common_utils/parameter.h"

struct SimpliestTreeStructureParameters : public ParametersSet
{
    int max_depth = 4;
    std::vector<float> branch_len = {60,20,12,7};
    std::vector<float> branch_r = {1.75,0.75,0.2, 0.075};
    std::vector<float> branch_angle = {PI/6, PI/6, PI/4, PI/4};
    std::vector<int> branch_count = {20,10,8,10};
    std::vector<float> branching_start = {0.5,0.0,0,0};
    float leaves_count = 0.75;
    float leaf_size = 3;
    virtual ParametersSet *copy() override
    { 
        auto Ps = new SimpliestTreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
    virtual void get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                                    ParameterVariablesSet v_set = ParameterVariablesSet::ALL_VALUES) override
    {

    }
    virtual glm::vec3 get_tree_max_size() override
    {
        float len = 0;
        for (float &l : branch_len)
            len += l;
        return glm::vec3(len,len,len);
    }
    virtual void save_to_blk(Block &b) override;
    virtual void load_from_blk(Block &b) override;
    virtual void write_parameter_list(ParameterList &list) override;
    virtual void read_parameter_list(ParameterList &list) override;
};

class SimpliestTreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override {};
private:
    void create_tree(Tree *tree, glm::vec3 pos, SimpliestTreeStructureParameters &params);
    void create_branch(Tree *tree, Branch *branch, glm::vec3 start_pos, glm::vec3 base_dir, glm::vec3 normal, int level, 
                       SimpliestTreeStructureParameters &params);
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    virtual bool iteration_method_implemented() {return false;}
    
    std::vector<glm::vec3> tree_positions;
    std::vector<TreeTypeData *> types;
};