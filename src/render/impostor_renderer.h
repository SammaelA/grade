#pragma once
#include "tree_utils/impostor.h"
#include "grove_renderer.h"

class ImpostorRenderer
{
    public:
    friend class GroveRenderer;
    ImpostorRenderer(ImpostorsData *data = nullptr);
    ~ImpostorRenderer();
    void render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light, 
                glm::mat4 &shadow_tr, GLuint shadow_tex,
                glm::vec3 camera_pos = glm::vec3(0,0,0),
                glm::vec4 screen_size = glm::vec4(800,600,1/800,1/600), 
                bool to_shadow = false,
                GroveRendererDebugParams dbgpar = GroveRendererDebugParams());
    private:
    GLuint slicesBuffer = 0;
    GLuint impostorsDataBuffer = 0;
    ImpostorsData *data= nullptr;

    Shader impostorRenderer;
    Shader impostorRendererInstancing;
    
    std::vector<Model *> models;
    std::vector<uint> offsets;

    float hth = 0,vth = 0;
};