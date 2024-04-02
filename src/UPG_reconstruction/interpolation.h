#pragma once 

#pragma once

#include <vector>
#include <iostream>
#include <random>
#include <omp.h>
#include "LiteMath.h"
#include <time.h>
#include <cmath>
#include "precompute_mat.h"


namespace interpolation
{
    const std::vector<float> create_A(const std::vector<LiteMath::float3>& X);
    void QR(const std::vector<float>& M, const size_t &size, std::vector<float> &Q, std::vector<float> &R);
    float matrix_norm(const std::vector<float>& A1, const std::vector<float>& A2);
    std::vector<float> calc_qr_coefs(const std::vector<float> &Q, const std::vector<float> &R, const std::vector<float> &b);
    std::vector<float> mul_qr(const std::vector<float> &Q, const std::vector<float> &R, size_t size);
    float perform_interpolation(const std::vector<float> &coefs, const LiteMath::float3 &pos);
    void householder_qr(const std::vector<float>& M, const size_t &size, std::vector<float> &Q, std::vector<float> &R);
    std::vector<float> calc_coefs(const std::vector<float> &b);
    float calc_interpolation(const std::vector<float> &coefs, const LiteMath::float3 &Point);

    const int interpolation_power = 4;
};