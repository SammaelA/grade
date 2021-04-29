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

class TerrainRenderer
{
public:
    TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step);
    ~TerrainRenderer();
    void render(glm::mat4 prc, glm::mat4 shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                                 DirectedLight &light);
private:
    float base_height = 0.0;
    Model *flat_terrain;
    Shader terrain;
};