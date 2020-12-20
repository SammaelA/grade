#pragma once
#include <list>
#include <vector>
#include "volumetric_occlusion.h"
#include "parameter.h"
class Texture;
class BillboardCloud;
class Shader;
struct Parameters
{
    glm::vec3 pos = glm::vec3(0,0,0);
    int max_depth = 4;
    int max_segments = 20;
    int max_branching = 2;
    int iterations = 60;

    float min_feed = 10;
    float max_feed = 10;
    float feed_v_mult = 0.01;

    float dir_conserv = 1.0; 
    float spread = 4;
    float phototrop = 1.0;
    float gravitrop = 0.8;

    float seg_dir_conserv = 1.0; 
    float seg_spread = 0.0;
    float seg_phototrop = 0.075;
    float seg_gravitrop = 0.1;
    float seg_bend = 0.075;
    float seg_bend_pow = 2.0;

    float base_r = 1.0;
    float r_level_mult = 0.4;
    float r_split_save_pow = 2.7;
    float r_feed_pow = 0.33;

    float scale = 3;
    float seg_len_mult = 5;
    float leaf_size_mult = 1.25;
    float len_level_mult = 0.45;

    float split_rate = 0.5; 
    float seg_split_rate = 0.1;
    float base_branch_feed = 1000;
    float base_seg_feed = 200;

    float light_precision = 1;
    float branch_grow_decrease_q = 0.9;
    float segment_grow_decrease_q = 0.09;

    float min_branching_chance = 0.0;
    float max_branching_chance = 1.0;
    float branching_power = 0.5;
    
    float feed_distribution_min_weight = 0.07;
    float feed_distribution_d_weight = 0.1;
    float top_growth_bonus = 0.15;

    float branch_removal = 0.75;
};
struct Joint;
struct Segment;
struct Branch;
struct Leaf;
struct Segment
{
    Joint *next;
    glm::vec3 begin;
    glm::vec3 end;
    float rel_r_begin;
    float rel_r_end;

};
struct Joint
{
    enum JointType {NONE, END, MIDDLE, FORK,LEAF};
    JointType type;
    Leaf *leaf = nullptr;
    glm::vec3 pos;
    std::list<Branch*> childBranches;
    int max_branching = 0;
    float light;
};
struct Branch
{
    int level;
    int base_seg_n;
    int max_seg_count;
    std::list<Segment> segments;
    std::list<Joint> joints;
    float light;
    float size;
    float base_r;
    bool dead = false;
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
    Leaf *new_leaf() {leaves.push_back(Leaf()); return &leaves.back();}
    void clear_removed();
};
struct BranchHeap
{
    std::list<Branch> branches;
    Branch *new_branch() {branches.push_back(Branch()); return &branches.back();}
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
    std::vector<BranchHeap*> branchHeaps;
    LeafHeap *leaves;
    glm::vec3 pos;
    TreeStructureParameters params;
    Branch *root;
    int iter = 0;
    uint id = 0;
    LightVoxelsCube *voxels;
    Texture *wood = nullptr;
    Texture *leaf = nullptr;
    std::vector<BillboardCloud *> billboardClouds;
    std::vector<Model *> models;
    void render(Shader &defaultShader, int cloudnum, glm::mat4 projcam);
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