#pragma once
#include "common_utils/field_2d.h"
#include "common_utils/utility.h"
#include "tinyEngine/texture.h"
class DirectedLight;
class Heightmap : public Field_2d
{
public:
    Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size);
    Heightmap(glm::vec3 pos, int w, int h);
    Heightmap() : Field_2d() {};
    float get_height(glm::vec3 pos);
    float get_height_simple(glm::vec3 pos);
    void random_generate(float base, float min, float max);
    void load_from_image(float base, float min, float max, std::string texture_name);
    glm::vec2 get_height_range() {return get_range();}
    glm::vec2 get_grad(glm::vec3 pos);

};
class HeightmapTex 
{
public:
    Texture &get() {return hmtex;}
    HeightmapTex(Heightmap &heightmap, int w = 1024, int h = 1024);
    ~HeightmapTex();
private:
    Texture hmtex;
    float base_value = 0;
    float size_x,size_y,height;
};