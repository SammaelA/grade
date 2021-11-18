#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
struct Scene
{
    Scene() {};
    ~Scene() {if (heightmap) delete heightmap;}
    Heightmap *heightmap = nullptr;
    GrovePacked grove;
};