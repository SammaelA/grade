#pragma once
#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"
#include <unordered_set>
#include "volumetric_occlusion.h"
#include "tinyEngine/utility/model.h"
#include "billboard_cloud.h"
#include "tree.h"
class DebugVisualizer;
class TreeGenerator
{
    public:
    TreeGenerator(Tree &t):curTree(t), curParams() {};
    bool random_tree(Tree &t);
    void plant_tree(Tree &t, TreeStructureParameters params);
    void grow_tree(Tree &t);
    void create_tree(Tree &t, TreeStructureParameters params, DebugVisualizer &debug);
    Tree &curTree;
    Branch *root;
    Branch *test;
    TreeStructureParameters curParams;
    LightVoxelsCube *voxels;
    bool random_iterate(Branch *b);
    void grow_branch(Branch *b, float feed);
    void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void seg_vertexes_to_model(SegmentVertexes &sv, Model *m);
    bool tree_to_model(Tree &t, bool leaves, bool debug = false);
    void leaf_to_model(Leaf &l, Model *m);
    void recursive_branch_to_model(Branch &b, Model *m, bool leaves);
    void branch_to_model(Branch &b, Model *m, bool leaves);
    void segment_to_model(Segment &s, Model *m, bool leaves);
    void joint_to_model(Joint &j, Model *m, bool leaves);
    void new_joint(Branch *b,Joint &j);
    void try_new_branch(Branch *base, Joint &j, Segment &s, bool from_end);
    void new_branch(Branch *base, Joint &j, Segment &s, glm::vec3 &M, bool from_end);
    void try_new_segment(Branch *base);
    void new_segment(Branch *base, glm::vec3 &M);
    void new_segment2(Branch *base, glm::vec3 &dir, glm::vec3 &pos);
    void set_seg_r(Branch *base, Segment &s);
    float calc_light(Joint& j);
    float calc_light(Branch *b);
    float calc_size(Joint& j);
    float calc_size(Branch *b);
    void distribute_feed(Branch *b);
    void remove_branch(Branch *b);
    void recalculate_thickness(Branch *b);
    LightVoxelsCube *create_light_voxels_cube(TreeStructureParameters params, glm::vec3 pos);
    glm::vec3 rand_dir();

    bool is_branch_productive(Branch *b);
};
