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

class ModelLoader;
class Visualizer
{
public:
friend class ModelLoader;
    Visualizer(Texture _tree_tex, Texture _leaves_tex, Shader *_tree_shader);
    Visualizer();
    void add_branch_layer(Tree &t, int layer, Mesh *m);
    void recursive_branch_to_model(Branch &b, Mesh *m, bool leaves, float scale = 1, int level_from = 0, 
                                   int level_to = 1000);
    void recursive_branch_to_model_fast(Branch &b, Mesh *m, bool leaves, float scale = 1, int level_from = 0, 
                                        int level_to = 1000);
    void leaf_to_model(Leaf &l, Mesh *m, float scale = 1);
    void packed_leaf_to_model(PackedLeaf &l, Mesh *m, glm::vec2 tc_zw = glm::vec2(1,0), bool need_tangent = false);
    void branch_to_model(Branch &b, Mesh *m, bool leaves);
    void packed_branch_to_model(PackedBranch &b, Mesh *m, bool leaves, float precision, glm::vec2 tc_zw = glm::vec2(1,0), bool need_tangents = false);
    void body_to_model(Body *b, Mesh *m, bool fixed_tc = false, glm::vec4 tc = glm::vec4(1,0,0,0));
    void heightmap_to_model(Heightmap &h, Mesh *m, glm::vec2 detailed_size, glm::vec2 full_size, float precision,
                            int LODs);
protected:
    void recursive_branch_to_model_fast_i(Branch &b, Mesh *m, bool leaves, float scale = 1, int level_from = 0, 
                                          int level_to = 1000);
    void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale);
    void get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, float rel_ring_pos,
                  std::vector<float> &mults, glm::vec3 p = glm::vec3(1,0,0));
    void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale);
    void seg_vertexes_to_model(SegmentVertexes &sv, Mesh *m, glm::vec2 tc_zw = glm::vec2(1,0), bool need_tangents = false);
    void segment_to_model(Segment &s, Mesh *m, bool leaves);
    void joint_to_model(Joint &j, Mesh *m, bool leaves);
    void box_to_model(Box *b, Mesh *m);
    void ellipsoid_to_model(Ellipsoid *b, Mesh *m, int sectors, int stacks, bool smooth = true);
    void cylinder_to_model(Cylinder *b, Mesh *m, int sectors);
    Texture tree_tex;
    Texture leaves_tex;
    Shader *tree_shader = nullptr;
};

class ModelLoader
{
public: 
    static void load_default_blk();
    Model *create_model_from_block(Block &bl, Texture &tex);
    Model *create_model_by_name(std::string name, Texture &tex);
    static std::string base_path;
private:
    Model *create_debug_box_model();
    Model *create_simple_grass_model();
    Model *load_model_from_obj(std::string name, Texture &tex);
};

extern Block obj_models_blk;