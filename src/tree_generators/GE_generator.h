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

class GETreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
    virtual bool iterate(LightVoxelsCube &voxels) override;
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    virtual bool iteration_method_implemented() override {return true;}
    static void set_joints_limit(int lim) {joints_limit = lim;}
private:
    static int joints_limit;
    static std::atomic<int> ids, t_ids;
    int iteration = 0;
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
    void set_occlusion_joint(Joint &j, float base_value, GETreeParameters &params, LightVoxelsCube &voxels);
};