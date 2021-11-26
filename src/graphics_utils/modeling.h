#pragma once

#include "core/tree.h"
#include "core/scene_models.h"
#include "core/grove.h"
#include "core/body.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/texture.h"
#include <vector>
#include "graphics_utils/volumetric_occlusion.h"
class Visualizer
{
public:
    Visualizer(Texture _tree_tex, Texture _leaves_tex, Shader *_tree_shader);
    Visualizer();
    void add_branch_layer(Tree &t, int layer, Model *m);
    void recursive_branch_to_model(Branch &b, Model *m, bool leaves, float scale = 1, int level_from = 0, int level_to = 1000);
    void leaf_to_model(Leaf &l, Model *m, float scale = 1);
    void packed_leaf_to_model(PackedLeaf &l, Model *m, glm::vec2 tc_zw = glm::vec2(1,0));
    void branch_to_model(Branch &b, Model *m, bool leaves);
    void packed_branch_to_model(PackedBranch &b, Model *m, bool leaves, int max_level, glm::vec2 tc_zw = glm::vec2(1,0));
    void body_to_model(Body *b, Model *m, bool fixed_tc = false, glm::vec4 tc = glm::vec4(1,0,0,0));
    void heightmap_to_model(Heightmap &h, Model *m, glm::vec2 detailed_size, glm::vec2 full_size, float precision,
                            int LODs);
protected:
    void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale);
    void get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, float rel_ring_pos,
                  std::vector<float> &mults, glm::vec3 p = glm::vec3(1,0,0));
    void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale);
    void seg_vertexes_to_model(SegmentVertexes &sv, Model *m, glm::vec2 tc_zw = glm::vec2(1,0));
    void segment_to_model(Segment &s, Model *m, bool leaves);
    void joint_to_model(Joint &j, Model *m, bool leaves);
    void box_to_model(Box *b, Model *m);
    void ellipsoid_to_model(Ellipsoid *b, Model *m, int sectors, int stacks, bool smooth = true);
    void cylinder_to_model(Cylinder *b, Model *m, int sectors);
    Texture tree_tex;
    Texture leaves_tex;
    Shader *tree_shader = nullptr;
};