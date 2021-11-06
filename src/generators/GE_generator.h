#pragma once

#include "../volumetric_occlusion.h"
#include "abstract_generator.h"
#include "../parameter.h"
#include "../volumetric_occlusion.h"
#include "../tinyEngine/utility/octree.h"
#include <vector>
#include <list>
struct GETreeParameters : public ParametersSet
{
    float lambda = 0.52;
    float k = 0.75;//part of joint that can create child branches
    int tau = 6;
    float ro = 1.0;
    float X0 = 2;
    float Xm = 100;
    float r = 0.4;
    int alpha = 4;
    float sigma = 0.5;
    float mu = 1.5;
    float nu = 1.0;
    float b_min = 1.8;
    float b_max = 2.2;
    float r_s = 0.1;

    float base_r = 0.025;
    int max_branches = 1;
    int occlusion_pyramid_d = 10;
    float r_pow = 2.2;
    int sp_points_base = 16;
    float branching_angle_min = 0;
    float branching_angle_max = PI/3;
    int max_iterations = 100;
    float leaf_size_mult = 3.5;
    float leaves_cnt = 0.0;
    int max_joints_in_branch = 8;
    float resource_mult = 7.0;
    float top_res_mult_base = 0.5;
    float top_res_mult_level_decrease = 0.5;
    float leaves_max_r = 2;//if radius in node > leaves_max_r*base_r leaf will not be created on this node

    virtual glm::vec3 get_tree_max_size() override
    {
        set_state(0);
        return ro*glm::vec3(2*Xm, 3.5*Xm, 2*Xm);
    }
};

class GETreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
    virtual bool iterate(LightVoxelsCube &voxels) override;
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    virtual bool iteration_method_implemented() override {return true;}
private:
    static int iteration, ids, t_ids;
    static GETreeParameters defaultParameters;
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
        int total_joints = 0;
        float base_r;
        Branch(){};
        Branch(int _level, glm::vec3 start_pos)
        {
            level = _level;
            joints = {Joint(start_pos,0,false)};
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
        Joint(glm::vec3 _pos, float _r, bool ch_b = true) 
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
    };

    struct SpaceColonizationData
    {
    public:
        void add(glm::vec3 pos)
        {
            positions.push_back(pos);
        }
        void prepare(LightVoxelsCube &voxels)
        {
            octree.create(voxels.get_bbox());
            octree.insert_vector(positions);
            positions.clear();
        }
        bool find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos, glm::vec3 dir, float angle,
                           glm::vec3 &best_pos, float &best_occ);
        void remove_close(glm::vec3 pos, float r);

        std::vector<glm::vec3> positions;
        Octree octree;
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
        Branch *base_branch;
        GrowthType gType;
        glm::vec3 prev_dir;
        float resource_left;
        GrowPoint(Joint *_j, Branch *_b, GrowthType _gt, glm::vec3 pd, float res_left)
        {
            joint = _j;
            base_branch = _b;
            gType = _gt;
            prev_dir = pd;
            resource_left = res_left;
        }
    };

    void create_tree_internal(Tree &t, GETreeParameters &params);
    void create_initial_trunk(Tree &t, GETreeParameters &params);
    void set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level);
    void convert(Tree &src, ::Tree &dst);
    void convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst);

    void calc_light(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params);
    void distribute_resource(Branch &b, GETreeParameters &params);
    void prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params, 
                                              std::vector<GrowPoint> &growth_points,
                                              SpaceColonizationData &sp_data,
                                              int max_growth_per_node);
    void grow_nodes(Tree &t, GETreeParameters &params, 
                    std::vector<GrowPoint> &growth_points,
                    SpaceColonizationData &sp_data,
                    LightVoxelsCube &voxels,
                    int max_growth_per_node);
    void remove_branches(Tree &t, Branch &b, GETreeParameters &params, LightVoxelsCube &voxels);
    void recalculate_radii(Tree &t, Branch &b, GETreeParameters &params);
    void add_SPCol_points_solid_angle(glm::vec3 pos, glm::vec3 dir, float r_max, int cnt, float min_psi, 
                                      SpaceColonizationData &sp_data);
    void set_occlusion(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params, float mul);
    void create_leaves(Branch &b, GETreeParameters &params, int level_from, LightVoxelsCube &voxels);
};