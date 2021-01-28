#pragma once
#include <list>
#include <vector>
#include "volumetric_occlusion.h"
#include "parameter.h"
#include "marks.h"
#include "grove.h"

class Texture;
class BillboardCloudRaw;
class Shader;
class Model;
struct Joint;
struct Segment;
struct Branch;
struct Leaf;
struct BranchHeap;
struct LeafHeap;

struct Segment
{
    Joint *next;
    glm::vec3 begin;
    glm::vec3 end;
    float rel_r_begin;
    float rel_r_end;
    Mark *mark = nullptr;
};
struct Joint
{
    enum JointType
    {
        NONE,
        END,
        MIDDLE,
        FORK,
        LEAF
    };
    JointType type;
    Leaf *leaf = nullptr;
    glm::vec3 pos;
    std::list<Branch *> childBranches;
    int max_branching = 0;
    float light;
    Mark *mark = nullptr;
};
struct Branch
{
    int id = 0;
    int level;
    int base_seg_n;
    int max_seg_count;
    std::list<Segment> segments;
    std::list<Joint> joints;
    float light;
    float size;
    float base_r;
    bool dead = false;
    glm::vec4 plane_coef;//plane ax+by+cz+d = 0 len(a,b,c) = 1
    void deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void transform(glm::mat4 &trans_matrix);
    void pack(PackedBranch &branch);
    Mark *mark = nullptr;
};
struct Leaf
{
    glm::vec3 pos;
    std::vector<glm::vec3> edges;
    bool dead = false;
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
struct TreeStructure
{
};
struct TreeData
{
};
struct Tree
{
    std::vector<BranchHeap *> branchHeaps;
    LeafHeap *leaves = nullptr;
    glm::vec3 pos;
    TreeStructureParameters params;
    Branch *root= nullptr;
    int iter = 0;
    uint id = 0;
    LightVoxelsCube *voxels= nullptr;
    Texture wood;
    Texture leaf;
    std::vector<BillboardCloudRaw *> billboardClouds;
    std::vector<Model *> models;
    void render(Shader &defaultShader, int cloudnum, glm::mat4 projcam);
    Tree();
    ~Tree();
};
struct VertexData
{
    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec2 tex_coord;
};
struct SegmentVertexes
{
    float ringsize = 0;
    std::vector<VertexData> smallRing;
    std::vector<VertexData> bigRing;
    Segment s;
};