#pragma once
#include <map>
#include <string>
#include <glm/glm.hpp>
#include <vector>
#include <list>

class Texture;
struct Joint;
struct Segment;
struct Branch;
struct Leaf;
struct BranchHeap;
struct LeafHeap;
struct Branch;
struct PackedBranch;
class Body;

struct Segment
{
    glm::vec3 begin;
    glm::vec3 end;
    float rel_r_begin;
    float rel_r_end;
    std::vector<float> mults;
};
struct Joint
{
    Leaf *leaf = nullptr;
    glm::vec3 pos;
    std::list<Branch *> childBranches;
};
struct Branch
{
    ushort type_id = 0;
    short level;
    std::list<Segment> segments;
    std::list<Joint> joints;
    glm::vec4 plane_coef;//plane ax+by+cz+d = 0 len(a,b,c) = 1
    glm::vec3 center_par;
    glm::vec3 center_self;
    
    int id = 0;
    int mark_a = 0;
    float mark_b = 0;

    void deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void norecursive_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void transform(glm::mat4 &trans_matrix);
    void pack(PackedBranch &branch);
    static float get_r_mult(float phi, std::vector<float> &mults);
};
struct Leaf
{
    glm::vec3 pos;
    std::vector<glm::vec3> edges;
    ushort type;
};
struct LeafHeap
{
    std::list<Leaf> leaves;
    Leaf *new_leaf()
    {
        leaves.push_back(Leaf());
        return &leaves.back();
    }
    void clear_removed();
    ~LeafHeap()
    {
        leaves.clear();
    }
    LeafHeap(){};

};
struct BranchHeap
{
    std::list<Branch> branches;
    Branch *new_branch()
    {
        branches.push_back(Branch());
        return &branches.back();
    }
    void clear_removed();
};

struct TreeTypeData
{
    TreeTypeData(int id, TreeStructureParameters params, std::string wood_tex_name, std::string leaf_tex_name);
    int type_id;
    TreeStructureParameters params;
    Texture wood;
    Texture leaf;
    int wood_id;
    int leaf_id;
    std::vector<Texture> additional_textures;

};
struct Tree
{
    std::vector<BranchHeap *> branchHeaps;
    LeafHeap *leaves = nullptr;
    glm::vec3 pos;
    ParameterSetWrapper params;
    Branch *root= nullptr;
    int iter = 0;
    uint id = 0;
    TreeTypeData *type = nullptr;
    LightVoxelsCube *voxels= nullptr;
    Texture wood;
    Texture leaf;
    std::vector<BillboardCloudRaw *> billboardClouds;
    std::vector<Model *> models;
    void render(Shader &defaultShader, int cloudnum, glm::mat4 projcam);
    Tree();
    ~Tree();
};
struct GroveGenerationData
{
    std::vector<TreeTypeData> types;
    int trees_count;
    int synts_count;
    int synts_precision;
    glm::vec3 pos;
    glm::vec3 size;
    std::string name;
    std::vector<Body *> obstacles;
};