#pragma once
#include "graphics_utils/terrain.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "generation/grove_generation_utils.h"

class Planter
{
public: 
    Planter(LightVoxelsCube *_voxels, Heightmap *_heightmap, GroveMask *_mask,
                  glm::vec3 center, glm::vec2 size,
                  float base_density, int max_saplings = 1000, float cell_size = 5);
    bool has_saplings() {return saplings_left > 0;}
    std::vector<glm::vec3> get_saplings();
private:
    LightVoxelsCube *voxels = nullptr;
    Heightmap *heightmap = nullptr;
    GroveMask *mask = nullptr;
    Field_2d density;
    Field_2d occlusion;

    int saplings_left = 0;
    glm::vec3 center;
    glm::vec2 size;
    float cell_size;
    std::vector<glm::vec3> saplings_planted;
};