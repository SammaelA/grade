#pragma once
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/shader.h"
#include "field_2d.h"
class DirectedLight;
class Heightmap : public Field_2d
{
public:
    Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size);
    Heightmap(glm::vec3 pos, int w, int h);

    float get_height(glm::vec3 pos);
    void random_generate(float base, float min, float max);
    glm::vec2 get_height_range() {return get_range();}
    glm::vec2 get_grad(glm::vec3 pos);

};
class HeightmapTex 
{
public:
    GLuint get() {return hmtex;}
    HeightmapTex(Heightmap &heightmap, int w = 1024, int h = 1024);
    ~HeightmapTex();
private:
    GLuint hmtex;
    float base_value = 0;
    float size_x,size_y,height;
};
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