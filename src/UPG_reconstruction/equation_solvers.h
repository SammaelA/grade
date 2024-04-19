#pragma once

#include <vector>
#include "LiteMath/LiteMath.h"

namespace solver
{
    std::vector<std::pair<float, float>> find_intervals(const std::vector<float>& coefs, const float& x1, const float& x2);
    float f(const std::vector<float>& coefs, const float& x);
}; // namespace solver
