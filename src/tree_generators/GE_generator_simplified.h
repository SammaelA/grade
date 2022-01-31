#pragma once

#include "graphics_utils/volumetric_occlusion.h"
#include "abstract_generator.h"
#include "common_utils/parameter.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "GE_generator_parameters.h"
#include <vector>
#include <list>
#include <atomic>


class GETreeGeneratorSimplified : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
    virtual bool iterate(LightVoxelsCube &voxels) override;
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    virtual bool iteration_method_implemented() override {return true;}
    virtual void set_seed(int _seed) override {gen = std::mt19937{_seed}; seed = _seed;}
    static void set_joints_limit(int lim) {joints_limit = lim;}
private:
    static int joints_limit;
    static std::atomic<int> ids, t_ids;
    static GETreeParameters defaultParameters;
    //std::random_device rd{};
    unsigned long seed = 0;
    std::mt19937 gen{0};
    std::uniform_real_distribution<double> d_ur = std::uniform_real_distribution<double>(0,1);
    int iteration = 0;
    struct Joint;
    struct Leaf; 
    struct Tree;

    std::vector<Tree> trees;

    struct Branch
    {
        std::list<Joint> joints;
        int level;
        bool alive = true;
        float total_resource = 0;
        float total_light = 0;
        int total_joints = 0;//how many joint this branch has totally (with subbranches)
        int distance_from_root = 0;//how many segnebts between first joint of this branch and root
        float base_r;
        glm::vec2 average_chb_dir = glm::vec2(0,0);
        Branch(){};
        Branch(int _level, glm::vec3 start_pos, int iteration)
        {
            level = _level;
            joints = {Joint(start_pos,0,iteration,false)};
        }
    };

    struct Leaf
    {
        glm::vec3 pos;
        std::vector<glm::vec3> edges;
    };

    struct Joint
    {
        int birth_time;
        glm::vec3 pos;
        float r;
        std::list<Branch> childBranches;
        Leaf leaf;//can be empty if edges.empty()
        float light = 0;
        float resource = 0;
        bool can_have_child_branches;
        Joint(glm::vec3 _pos, float _r, int iteration, bool ch_b) 
        {
            birth_time = iteration;
            pos = _pos; 
            r = _r; 
            childBranches = {};
            leaf = Leaf();
            resource = 0;
            light = 0;
            can_have_child_branches = ch_b;
        }
    };

    enum TreeStatus
    {
        SEED,
        GROWING,
        GROWN,
        DEAD
    };
    struct Tree
    {
        glm::vec3 pos;
        Branch root;
        int max_depth = 1;
        TreeTypeData *type = nullptr;
        TreeStatus status = SEED;
        int iteration = 0;
        int joints_total = 0;
    };
    
    enum GrowthType
    {
        END,
        BRANCHING,
        END_BRANCH,
        FINISHED
    };
    struct GrowPoint
    {
        Joint *joint;
        int joint_n;
        Branch *base_branch;
        GrowthType gType;
        glm::vec3 prev_dir;
        float resource_left;
        GrowPoint(Joint *_j, Branch *_b, GrowthType _gt, glm::vec3 pd, float res_left, int _joint_n)
        {
            joint = _j;
            base_branch = _b;
            gType = _gt;
            prev_dir = pd;
            resource_left = res_left;
            joint_n = _joint_n;
        }
    };

    void create_tree_internal(Tree &t, GETreeParameters &params);
    void create_initial_trunk(Tree &t, GETreeParameters &params);
    void set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level);
    void convert(Tree &src, ::Tree &dst);
    void convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst);

    void calc_light(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params);
    void distribute_resource(Branch &b, GETreeParameters &params, float res_mult);
    void prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params, 
                                              std::vector<GrowPoint> &growth_points,
                                              int max_growth_per_node);
    void grow_nodes(Tree &t, GETreeParameters &params, 
                    std::vector<GrowPoint> &growth_points,
                    LightVoxelsCube &voxels,
                    int max_growth_per_node);
    void remove_branches(Tree &t, Branch &b, GETreeParameters &params, LightVoxelsCube &voxels);
    void recalculate_radii(Tree &t, Branch &b, GETreeParameters &params);
    void set_occlusion(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params, float mul);
    void create_leaves(Branch &b, GETreeParameters &params, int level_from, LightVoxelsCube &voxels);
    void set_occlusion_joint(Joint &j, float base_value, GETreeParameters &params, LightVoxelsCube &voxels);
    bool find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos,
                       glm::vec3 dir, float angle,
                       glm::vec3 &best_pos, float &best_occ);
    inline float self_rand(double from = 0.0, double to = 1.0) 
    { 
        seed = seed * 1103515245 + 12345;
        float r =  ((unsigned int)(seed / 65536) % 32768) / 32768.0f;
        return from >= to ? from : from + r*(to - from);
    }
};