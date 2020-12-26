#pragma once
#include "tree.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/texture.h"
#include <vector>
class DebugVisualizer
{
public:
    DebugVisualizer(Texture *_tree_tex = nullptr, Shader *_tree_shader = nullptr);
    void render(glm::mat4 view_proj);
    void add_branch(Branch *b, glm::vec3 scale, glm::vec3 shift, int level = -1);
    void enable_all();
    void disable_all();
    ~DebugVisualizer() {};
private:
    Texture *tree_tex = nullptr;
    Shader *tree_shader = nullptr; 
    std::vector<Model *> debugModels;
    std::vector<int> currentModes;
};