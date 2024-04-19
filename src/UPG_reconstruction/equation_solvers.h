#pragma once

#include <vector>
#include "LiteMath/LiteMath.h"
#include <algorithm>
#include <iostream>


namespace solver
{
    void find_interval(const std::vector<float>& coefs, const float& x1, const float& x2, const int& n, float intervals[2]);
    float f(const std::vector<float>& coefs, const float& x, const int& n);
}; // namespace solver
