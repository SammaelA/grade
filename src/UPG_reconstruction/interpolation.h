#pragma once 

#include <vector>
#include <cstdlib>
#include <glm/glm.hpp>

namespace interpolation
{
    float bilinear(const float &tx, const float &ty, const float &c00, const float &c10, const float &c01, const float &c11);
    float trilinear(const std::vector<float>& sdf_model, glm::vec3& p, size_t dim_size);
};