#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
#include "grass.h"
struct TreeTypeData;
struct Scene
{
    Scene() {};
    ~Scene() {if (heightmap) delete heightmap;}
    Heightmap *heightmap = nullptr;
    GrovePacked grove;
    GrassPacked grass;
    std::vector<TreeTypeData> tree_types;
};