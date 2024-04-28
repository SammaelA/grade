#pragma once

#include <vector>
#include "LiteMath/LiteMath.h"
#include <algorithm>
#include <iostream>


namespace solver
{
    void find_interval(const float* coefs, const float& x1, const float& x2, const int& n, float intervals[2]);
    float f(const float* coefs, const float& x, const int& n);
    void fd(const std::vector<float>& coefs, const float& x, const float& last_x, const int& n, float res[2]);
    float find_root();
    void polinomMiltiplier(const float* p1, const float* p2, float* res, const int& n1, const int& n2, const float& polinom_coef);
    //  before call this func need to set Ray coords from world position to voxel space
    void coefsDecrease(const std::vector<float>& coefs, const LiteMath::float3& P, const LiteMath::float3& D, float* res);

    //  For test
    float calc_test_res(const std::vector<float>& coefs, const LiteMath::float3& P, const LiteMath::float3& D, const float& t);
};
