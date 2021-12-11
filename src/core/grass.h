#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "graphics_utils/texture_atlas.h"
#include "save_utils/blk.h"

struct GrassType
{
    int id = -1;

    Texture texture;
    float patch_size = 25;
    float patch_size_std_dev = 10;
    float patch_density = 1;
    float patch_density_std_dev = 0.5;
    float plant_size = 3;
    float plant_size_std_dev = 1.5;
    float light_sensivity = 1;
    float push = 1;
    std::string model_name = "simple_grass";
    void load_from_blk(Block &block);
    GrassType();
};

struct GrassInstanceData
{
    glm::vec3 pos;
    float size;
    float rot_y;
    GrassInstanceData(glm::vec3 _pos, float _size, float _rot_y)
    {
        pos = _pos;
        size = _size;
        rot_y = _rot_y;
    }
};

struct GrassPacked
{
    std::vector<GrassType> used_grass_types;
    TextureAtlas grass_textures;
    std::vector<std::pair<int, std::vector<GrassInstanceData>>> grass_instances;
    //each member is a type: int is texture id in atlas, second is list of instances for this type
};