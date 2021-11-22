#pragma once
#include "terrain_renderer.h"
#include "core/grass.h"

class GrassRenderer
{
public:
    GrassRenderer();
    void render(glm::mat4 &projection, glm::mat4 &view, glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
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
    GrassRenderer2(GrassPacked &data);
    ~GrassRenderer2();
    void render(glm::mat4 &projection, glm::mat4 &view, glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow = false);
private:
    Shader grass;
    Shader grassShadow;
    std::vector<Model *> models;
    std::vector<int> inst_offsets;
    std::vector<int> inst_counts;
    GLuint instances_buffer;
    TextureAtlas &grass_atlas;
};