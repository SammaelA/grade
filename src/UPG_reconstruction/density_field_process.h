#pragma once

#include "obj_utils.h"
#include <string>
#include <iostream>
#include <ctime>

#define samplePointsCount 10
#define density_dim 128

namespace df
{
    struct Voxel 
    {
        float size = 0;
        glm::vec3 samplingPoints[samplePointsCount];
        size_t inside_count = 0;
        float density = 0;
    };

    struct VoxelGrid 
    {
        Voxel *voxels;
        size_t dimension = density_dim;
        glm::vec3 bounds[2] = { glm::vec3(-5.f,-5.f,-5.f), glm::vec3(5.f, 5.f, 5.f) };
    };

    float* create_density_field(const std::vector<float>& model);
    void erase(float* density, float trashold = 0);
    void create_sdf(float* density);
    float* pipeline(const std::vector<float>& model);
};