#pragma once
#include "tree.h"
#include "grove.h"
#include "body.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/texture.h"
#include <vector>
#include "volumetric_occlusion.h"
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
protected:
    void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, float rel_ring_pos,
                  std::vector<float> &mults, glm::vec3 p = glm::vec3(1,0,0));
    void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
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
class DebugVisualizer : public Visualizer
{
public:
    static const int MAX_RENDER_MODE = 2;
    DebugVisualizer();
    DebugVisualizer(Texture _tree_tex, Shader *_tree_shader = nullptr);
    void render(glm::mat4 view_proj, int mode = -1);
    void add_branch_debug(Branch *b, glm::vec3 scale, glm::vec3 shift, int level = -1);
    void enable_all();
    void disable_all();
    void add_bodies(Body *b_ptr, int count);
    void visualize_light_voxels(LightVoxelsCube *voxels, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1));
    void visualize_light_voxels(LightVoxelsCube *voxels,glm::vec3 pos, glm::vec3 size, glm::vec3 step, float dot_size, 
                                float threshold = 0, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1),
                                int mip = 0);
    ~DebugVisualizer();
    std::vector<Model *> debugModels;
    Shader debugShader;
    Shader bodyShader; 
    
    DebugVisualizer& operator=(const DebugVisualizer& dv);

private:
    void branch_to_model_debug(Branch *b, int level, Model &m);
    std::vector<int> currentModes;
};