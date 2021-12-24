#pragma once
#include <list>
#include <vector>
#include <string>
#include "common_utils/utility.h"
#include "../tinyEngine/texture.h"
#include "common_utils/parameter.h"

class Texture;
struct Joint;
struct Segment;
struct Branch;
struct Leaf;
struct BranchHeap;
struct LeafHeap;
struct PackedBranch;

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
    short mark_A;
};
struct Branch
{
    bool dead = false;
    int self_id = 0;
    int id = 0;
    ushort type_id = 0;
    short level;
    std::list<Segment> segments;
    std::list<Joint> joints;
    glm::vec4 plane_coef;//plane ax+by+cz+d = 0 len(a,b,c) = 1
    glm::vec3 center_par;
    glm::vec3 center_self;
    int mark_A;
    int mark_B;
    void deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void norecursive_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void transform(glm::mat4 &trans_matrix, float r_transform = 1.0);
    void pack(PackedBranch &branch);
    void mark_dead();
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
    LeafHeap& operator=(LeafHeap&&h)
    {
        leaves = std::move(h.leaves);
    }

};
struct BranchHeap
{
    std::list<Branch> branches;
    Branch *new_branch()
    {
        branches.push_back(Branch());
        return &branches.back();
    }
    BranchHeap(){};
    BranchHeap& operator=(BranchHeap&&h)
    {
        branches = std::move(h.branches);
    }
    void clear_removed();
};

struct TreeTypeData
{
    TreeTypeData();
    TreeTypeData(int id, ParametersSet *params, std::string wood_tex_name, std::string leaf_tex_name);
    ~TreeTypeData();
    int type_id;
    ParametersSet *params = nullptr;
    Texture wood;
    Texture leaf;
    int wood_id;
    int leaf_id;
    std::vector<Texture> additional_textures;
    std::string generator_name = "default";

};
struct Tree
{
    std::vector<BranchHeap *> branchHeaps;
    LeafHeap *leaves = nullptr;
    glm::vec3 pos = glm::vec3(0,0,0);
    Branch *root= nullptr;
    uint id = 0;
    TreeTypeData *type = nullptr;
    bool valid = false;
    Tree() {};
    ~Tree() {clear();};
    void clear()
    {
        if (leaves)
        {
            delete leaves;
            leaves = nullptr;
        }
        for (int i=0;i<branchHeaps.size();i++)
        {
            delete branchHeaps[i];
            branchHeaps[i] = nullptr;
        }
        branchHeaps.clear();
        root = nullptr;
        type = nullptr;
    }
};

enum Quality
{
    LOW_AS_F = 128,
    ULTRALOW = 256,
    LOW = 400,
    MEDIUM = 600,
    HIGH = 800,
    ULTRA = 1024
};
