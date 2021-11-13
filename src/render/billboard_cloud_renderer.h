#pragma once
#include "graphics_utils/billboard_cloud.h"
#include "grove_renderer.h"

struct BillboardCloudRenderer : Countable
{   
    friend class GroveRenderer;
    enum RenderMode
    {
        NOTHING,
        ONLY_SINGLE,
        ONLY_INSTANCES,
        BOTH
    };
    BillboardCloudRenderer(BillboardCloudData *data = nullptr);
    ~BillboardCloudRenderer();
    void render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light,
                glm::mat4 &shadow_tr, GLuint shadow_tex,
                glm::vec3 camera_pos = glm::vec3(0,0,0),
                glm::vec4 screen_size = glm::vec4(800,600,1/800,1/600),
                bool to_shadow = false,
                GroveRendererDebugParams dbgpar = GroveRendererDebugParams());
    void set_render_mode(RenderMode m)
    {
        renderMode = m;
    }
    private:
    BillboardCloudData *data= nullptr;
    Shader rendererToTexture;
    Shader billboardRenderer;
    Shader billboardRendererInstancing;
    Model *cloud = nullptr;
    std::vector<Model *> instances;
    RenderMode renderMode = ONLY_SINGLE;

};