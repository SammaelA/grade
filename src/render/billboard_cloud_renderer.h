#pragma once
#include "tree_utils/billboard_cloud.h"
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
    BillboardCloudRenderer(const BillboardCloudData *data = nullptr);
    ~BillboardCloudRenderer();
    void render(MultiDrawRendDesc &mdrd, float4x4 &projection, float4x4 &view, DirectedLight &light,
                float4x4 &shadow_tr, GLuint shadow_tex,
                float3 camera_pos = float3(0,0,0),
                float4 screen_size = float4(800,600,1/800,1/600),
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