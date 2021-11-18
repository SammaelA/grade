#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
struct TreeTypeData;
struct Scene
{
    Scene() {};
    ~Scene() {if (heightmap) delete heightmap;}
    Heightmap *heightmap = nullptr;
    GrovePacked grove;
    std::vector<TreeTypeData> tree_types;
};