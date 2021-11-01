#pragma once

#include "../volumetric_occlusion.h"
#include "abstract_generator.h"
#include "../parameter.h"
#include "../volumetric_occlusion.h"
#include "../tinyEngine/utility/octree.h"
#include <vector>
#include <list>
struct GETreeParameters /*: public ParametersSet*/
{
    float lambda = 0.52;
    float k = 0.5;//part of joint that can create child branches
    int tau = 6;
    float ro = 1.5;
    float X0 = 2;
    float Xm = 40;
    float r = 0.34;
    int alpha = 4;
    float sigma = 0.5;
    float mu = 1.5;
    float nu = 0.5;
    float b_min = 1.8;
    float b_max = 2.2;
    float r_s = 0.25;

    float base_r = 0.04;
    int max_branches = 2;
    int occlusion_pyramid_d = 10;
    float r_pow = 2.5;
    int sp_points_base = 10;
    float branching_angle = PI/4;
    int max_iterations = 20;
};

class GETreeGenerator : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
private:
    static int iteration, ids, t_ids;
    
    struct Joint;
    struct Leaf; 
    struct Branch
    {
        std::list<Joint> joints;
        int level;
        bool alive = true;
        float total_resource;
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
        float resource;
        bool can_have_child_branches;
        Joint(glm::vec3 _pos, float _r, bool ch_b = true) 
        {
            birth_time = iteration;
            pos = _pos; 
            r = _r; 
            childBranches = {};
            leaf = Leaf();
            resource = 0;
            can_have_child_branches = ch_b;
        }
    };

    struct Tree
    {
        glm::vec3 pos;
        Branch root;
        int max_depth = 1;
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
        /*
        {
            float r_sqr = r*r;
            best_occ = 1000;
            float cs = cos(angle);
            for (auto &p : positions)
            {
                if (dot(p-pos,p-pos) <= r_sqr && dot(normalize(p-pos),dir) > cs)
                {
                    float occ = voxels.get_occlusion(p);
                    //logerr("%f occ", occ);
                    if (occ < best_occ)
                    {
                        best_occ = occ;
                        best_pos = p;
                    }
                }
            }
            if (best_occ >= 1000)
                return false;
            else 
                return true;
        }
        */
        void remove_close(glm::vec3 pos, float r);
        /*
        {
            float r_sq = r*r;
            auto it = positions.begin();
            while (it != positions.end())
            {
                float d = glm::dot(pos - *it, pos - *it);
                if (d < r_sq)
                {
                    it = positions.erase(it);
                }
                else
                    it++;
            }
        }*/

        std::vector<glm::vec3> positions;
        Octree octree;
    };
    
    enum GrowthType
    {
        END,
        BRANCHING,
        FINISHED
    };
    struct GrowPoint
    {
        Joint *joint;
        Branch *base_branch;
        GrowthType gType;
        glm::vec3 prev_dir;
        GrowPoint(Joint *_j, Branch *_b, GrowthType _gt, glm::vec3 pd)
        {
            joint = _j;
            base_branch = _b;
            gType = _gt;
            prev_dir = pd;
        }
    };

    void create_tree_internal(Tree &t, GETreeParameters &params);
    void create_initial_trunk(Tree &t, GETreeParameters &params);
    void set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level);
    void convert(Tree &src, ::Tree &dst, GroveGenerationData &ggd);
    void convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst);

    void calc_light(Tree &t, Branch &b, LightVoxelsCube &voxels);
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
    void set_occlusion(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params);
};