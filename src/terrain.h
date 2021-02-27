#pragma once
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/shader.h"
class Heightmap
{
public:
    Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size);
    Heightmap(glm::vec3 pos, int w, int h);
    ~Heightmap();
    float get_height(glm::vec3 pos);
    void random_generate(float base, float min, float max);
private:
    float get(int x, int y);
    int w, h;
    glm::vec3 pos;
    glm::vec2 size;
    float cell_size = 1;
    float *data;
    float base_height;

};

class TerrainRenderer
{
public:
    TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step);
    ~TerrainRenderer();
    void render(glm::mat4 prc);
private:
    float base_height = 0.0;
    Model *flat_terrain;
    Shader terrain;
};