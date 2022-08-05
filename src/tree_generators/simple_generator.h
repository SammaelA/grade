#pragma once
#include "abstract_generator.h"
#include "common_utils/parameter.h"
#include "mygen_parameters.h"
struct SimpleTreeStructureParameters : public ParametersSet
{
    int max_depth = 4;
    glm::vec4 segment_size = glm::vec4(3.0,1.5,0.5,0.5);
    glm::vec4 segment_count = glm::vec4(15, 12, 10, 15);
    glm::vec4 base_thickness = glm::vec4(2,1.4,0.6,0.4);
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

    virtual void save_to_blk(Block &b) override
    {
      ParameterList list;
      save_load_define(SaveLoadMode::BLK_SAVE, b, list);
    }
    virtual void load_from_blk(Block &b) override
    {
      ParameterList list;
      save_load_define(SaveLoadMode::BLK_LOAD, b, list);
    }
    virtual void RW_parameter_list(bool write, ParameterList &list) override 
    {
      Block b;
      save_load_define(write ? SaveLoadMode::PAR_LIST_LOAD : SaveLoadMode::PAR_LIST_SAVE, b, list);
    }
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
    virtual glm::vec3 get_tree_max_size() override
    {
      return 2.0f*glm::vec3(segment_size[0]*segment_count[0]);
    }
    virtual ParametersSet *copy() override
    { 
        auto Ps = new SimpleTreeStructureParameters();
        *Ps = *this;
        return Ps;
    };
};

class SimpleTreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
private:
    void create_tree(Tree *tree, glm::vec3 pos);
    void create_branch(Tree *tree, Branch *branch, glm::vec3 start_pos, glm::vec3 base_dir, glm::vec3 normal, int level, 
                       float base_r, float leaves_chance);
    SimpleTreeStructureParameters params;
    GroveGenerationData ggd;
    Heightmap *h;
    std::vector<glm::vec3> tree_positions;
    std::vector<TreeTypeData *> types;
};