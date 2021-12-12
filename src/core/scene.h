#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
#include "grass.h"
#include "common_utils/bvh.h"

struct TreeTypeData;
struct Scene
{
    enum ObjCategories
    {
        UNKNOWN,
        ERROR,
        SIMPLE_OBJECT
    };
    struct InstancedModel
    {
        std::string name;
        Texture tex;
        Model *model;
        std::vector<glm::mat4> instances;
        InstancedModel();
    };
    Scene(){};
    ~Scene();
    Heightmap *heightmap = nullptr;
    std::vector<InstancedModel> instanced_models;

    GrovePacked grove;
    GrassPacked grass;
};