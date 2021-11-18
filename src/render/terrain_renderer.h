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
                DirectedLight &light, bool to_shadow = false);
    Model *flat_terrain;
private:
    float base_height = 0.0;
    Texture terrain_tex;
    Shader terrain;
    Shader terrainShadow;
};