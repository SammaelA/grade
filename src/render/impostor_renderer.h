#pragma once
#include "tree_utils/impostor.h"
#include "grove_renderer.h"

class ImpostorRenderer
{
    public:
    friend class GroveRenderer;
    ImpostorRenderer(ImpostorsData *data = nullptr);
    ~ImpostorRenderer();
    void render(MultiDrawRendDesc &mdrd, float4x4 &projection, float4x4 &view, DirectedLight &light, 
                float4x4 &shadow_tr, GLuint shadow_tex,
                float3 camera_pos = float3(0,0,0),
                float4 screen_size = float4(800,600,1/800,1/600), 
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