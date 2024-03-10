#pragma once
#include "graphics_utils/terrain.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
class TerrainRenderer
{
public:
    TerrainRenderer(Heightmap &h, float3 pos, float2 size, float2 step);
    ~TerrainRenderer();
    void render(float4x4 projection, float4x4 view, float4x4 shadow_tr, GLuint shadow_tex, float3 camera_pos,
                DirectedLight &light, bool to_shadow = false, 
                int debug_type = 0, float4 grid_params = float4(0,0,1,1), float4 debug_tex_scale = float4(0,0,0.01, 0.01),
                Texture *debug_tex_ptr = nullptr);
    Model *flat_terrain;
private:
    Texture terrain_tex1;
    Texture terrain_tex2;
    Texture terrain_tex3;
    Texture perlin;
    Shader terrain;
    Shader terrainShadow;
};