#pragma once

#include "../volumetric_occlusion.h"
#include "abstract_generator.h"
#include "../parameter.h"

struct GETreeParameters /*: public ParametersSet*/
{
    float lambda = 0.52;
    float k = 0.5;
    int tau = 6;
    float ro = 0.1;
    float X0 = 2;
    float Xm = 100;
    float r = 0.34;
    int alpha = 4;
    float sigma = 0.5;
    float mu = 1.5;
    float nu = 0.5;
    float b_min = 1.8;
    float b_max = 2.2;
    float r_s = 0.25;

    float base_r = 0.05;
    float seg_len = 1;
    int max_branches = 1;
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
        std::vector<Joint> joints;
        int level;
        bool alive = true;
        float total_resource;
        int total_joints = 0;
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
        std::vector<Branch> childBranches;
        Leaf leaf;//can be empty if edges.empty()
        float resource;

        Joint(glm::vec3 _pos, float _r) 
        {
            birth_time = iteration;
            pos = _pos; 
            r = _r; 
            childBranches = {};
            leaf = Leaf();
            resource = 0;
        }
    };

    struct Tree
    {
        glm::vec3 pos;
        Branch root;
        int max_depth = 1;
    };

    void create_tree_internal(Tree &t, GETreeParameters &params);
    void create_initial_trunk(Tree &t, GETreeParameters &params);
    void set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level);
    void convert(Tree &src, ::Tree &dst, GroveGenerationData &ggd);
    void convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst);
};