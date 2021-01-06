#pragma once
#include "tree.h"
#include "grove.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/texture.h"
#include <vector>
class Visualizer
{
public:
    Visualizer(Texture *_tree_tex, Texture *_leaves_tex, Shader *_tree_shader);
    Visualizer(){};
    void add_branch_layer(Tree &t, int layer, Model *m);
    void recursive_branch_to_model(Branch &b, Model *m, bool leaves);
    void leaf_to_model(Leaf &l, Model *m);
    void branch_to_model(Branch &b, Model *m, bool leaves);
    void packed_branch_to_model(PackedBranch &b, Model *m, bool leaves);
    void set_params(TreeStructureParameters &params) { curParams = params; }

protected:
    void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos);
    void seg_vertexes_to_model(SegmentVertexes &sv, Model *m);
    void segment_to_model(Segment &s, Model *m, bool leaves);
    void joint_to_model(Joint &j, Model *m, bool leaves);
    TreeStructureParameters curParams;
    Texture *tree_tex = nullptr;
    Texture *leaves_tex = nullptr;
    Shader *tree_shader = nullptr;
};
class DebugVisualizer : public Visualizer
{
public:
    DebugVisualizer(Texture *_tree_tex = nullptr, Shader *_tree_shader = nullptr);
    void render(glm::mat4 view_proj);
    void add_branch_debug(Branch *b, glm::vec3 scale, glm::vec3 shift, int level = -1);
    void enable_all();
    void disable_all();
    ~DebugVisualizer();
    std::vector<Model *> debugModels;

private:
    void branch_to_model_debug(Branch *b, int level, Model &m);
    std::vector<int> currentModes;
};