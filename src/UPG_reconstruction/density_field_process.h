#pragma once

#include "obj_utils.h"
#include <string>
#include <iostream>
#include <ctime>
#include <vector>

namespace df
{
    //? Maybe this structure is not necessary
    struct Voxel 
    {
        std::vector<glm::vec3> samplingPoints;
        float density = 0;
    };

    struct VoxelGrid 
    {
        //? Maybe this array is not necessary
        //? std::vector<Voxel> voxels;
        size_t dimension;
        size_t samplePointsCount;
        glm::vec3 bounds[2];

        //?  Init Voxel Grid using dimension size, sample points count and grid bounds
        VoxelGrid(size_t dim, size_t samplePointsCount, glm::vec3 bounds[2])
        {
            this->dimension = dim;
            this->samplePointsCount = samplePointsCount;
            this->bounds[0] = bounds[0];
            this->bounds[1] = bounds[1];
        };
    };

    std::vector<float> create_density_field(const std::vector<float>& model, const VoxelGrid& grid);
    void erase(std::vector<float>& density, const VoxelGrid& grid, const float& trashold = 0);
    std::vector<float> create_sdf(const std::vector<float>& density, const VoxelGrid& grid);
    std::vector<float> pipeline(const std::vector<float>& model);
    glm::vec3 closest_point_triangle(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c);
};