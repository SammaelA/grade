#pragma once
#include "graphics_utils/terrain.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
class TerrainRenderer
{
public:
    TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step);
    ~TerrainRenderer();
    void render(glm::mat4 projection, glm::mat4 view, glm::mat4 shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                DirectedLight &light, bool to_shadow = false, 
                int debug_type = 0, glm::vec4 grid_params = glm::vec4(0,0,1,1), glm::vec4 debug_tex_scale = glm::vec4(0,0,0.01, 0.01),
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