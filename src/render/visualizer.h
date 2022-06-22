#pragma once
#include "graphics_utils/modeling.h"
#include <vector>
#include "common_utils/bbox.h"
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

namespace visualizer
{
    Model *visualize_light_voxels(LightVoxelsCube *voxels, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1));
    Model *visualize_light_voxels(LightVoxelsCube *voxels,glm::vec3 pos, glm::vec3 size, glm::vec3 step, float dot_size, 
                                  float threshold = 0, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1),
                                  int mip = 0);
    Model *visualize_aabb(AABB &box, glm::vec3 &color);
    Model *visualize_aabb(::std::vector<AABB> &boxes, ::std::vector<glm::vec3> &colors);
};