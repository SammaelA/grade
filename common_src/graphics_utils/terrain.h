#pragma once
#include "common_utils/field_2d.h"
#include "common_utils/utility.h"
#include "tinyEngine/texture.h"
struct DirectedLight;
class Heightmap : public Field_2d
{
public:
    Heightmap(float3 pos, float2 size, float cell_size);
    Heightmap(float3 pos, int w, int h);
    Heightmap() : Field_2d() {};
    float get_height(float3 pos);
    float get_height_simple(float3 pos);
    void random_generate(float base, float min, float max);
    void load_from_image(float base, float min, float max, std::string texture_name);
    float2 get_height_range() {return get_range();}
    float2 get_grad(float3 pos);

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
};