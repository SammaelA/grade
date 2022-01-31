#pragma once 
#include "bbox.h"
#include <list>
#include <functional>
#include <vector>
#include "graphics_utils/volumetric_occlusion.h"
#include <glm/glm.hpp>

class VoxList
{
public:
    static constexpr int base_points_cnt = 16;
    void create(LightVoxelsCube *voxels);
private:
    struct Vox
    {
        int points_cnt = 0;
        glm::vec3 points[base_points_cnt];  
        std::vector<glm::vec3> more_points; 
    };
};