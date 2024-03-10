#pragma once
#include "terrain_renderer.h"
#include "core/grass.h"

class GrassRenderer
{
public:
    GrassRenderer();
    void render(float4x4 &projection, float4x4 &view, float4x4 &shadow_tr, GLuint shadow_tex, float3 camera_pos,
                HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow = false);
private:
    Texture grass_base;
    Texture grass_tall;
    Texture perlin;
    Texture noise;
    Shader grass;
    Shader grassShadow;
    Model m;
};

class GrassRenderer2
{
public:
    GrassRenderer2(const GrassPacked &data);
    ~GrassRenderer2();
    void render(float4x4 &projection, float4x4 &view, float4x4 &shadow_tr, GLuint shadow_tex, float3 camera_pos,
                HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow = false);
private:
    Shader grass;
    Shader grassShadow;
    std::vector<Model *> models;
    std::vector<int> inst_offsets;
    std::vector<int> inst_counts;
    GLuint instances_buffer;
    const TextureAtlas &grass_atlas;
};