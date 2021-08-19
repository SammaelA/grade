#pragma once
#include <list>
#include <vector>
#include <string>
#include "volumetric_occlusion.h"
#include "parameter.h"
#include "mygen_parameters.h"
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
struct LeafHeap : Countable
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
struct BranchHeap : Countable
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
    TreeTypeData(int id, ParametersSet *params, std::string wood_tex_name, std::string leaf_tex_name);
    ~TreeTypeData();
    int type_id;
    ParametersSet *params = nullptr;
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
    ULTRALOW = 256,
    LOW = 400,
    MEDIUM = 600,
    HIGH = 800,
    ULTRA = 1024
};
enum GenerationTask
{
    GENERATE = 1,
    CLUSTERIZE = 2,
    BILLBOARDS = 4,
    IMPOSTORS = 8,
    IMPOSTOR_FULL_GROVE = 16,
    SYNTS = 32,
    MODELS = 64
};
#define ALL_GENERATION_TASKS (2*MODELS - 1)
#define MINIMUM_FOR_RENDER (GENERATE|MODELS|BILLBOARDS)
struct GroveGenerationData
{
    std::vector<TreeTypeData> types;
    int trees_count;
    int synts_count;
    int synts_precision;
    float clustering_max_individual_distance = 0.7;
    Quality impostor_quality = Quality::MEDIUM;
    Quality bill_1_quality = Quality::MEDIUM;
    Quality bill_2_quality = Quality::LOW;
    unsigned task = GENERATE | CLUSTERIZE | BILLBOARDS | IMPOSTORS | MODELS;
    glm::vec3 pos;
    glm::vec3 size;
    std::string name;
    std::vector<Body *> obstacles;
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