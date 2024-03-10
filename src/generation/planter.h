#pragma once
#include "graphics_utils/terrain.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "generation/grove_generation_utils.h"

class Planter
{
public: 
    Planter(LightVoxelsCube *_voxels, Heightmap *_heightmap, GroveMask *_mask, GroveMask *_biome_mask,
                  float3 center, float2 size,
                  float base_density, int max_saplings = 1000, float cell_size = 5);
    bool has_saplings() {return saplings_left > 0;}
    std::vector<float3> get_saplings();
private:
    LightVoxelsCube *voxels = nullptr;
    Heightmap *heightmap = nullptr;
    GroveMask *mask = nullptr;
    GroveMask *biome_mask = nullptr;
    Field_2d density;
    Field_2d occlusion;

    int saplings_left = 0;
    float3 center;
    float2 size;
    float cell_size;
    std::vector<float3> saplings_planted;
};