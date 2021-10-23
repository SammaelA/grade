#pragma once
#include <vector>

struct Hash
{
    std::vector<float> data;
    std::vector<int> start_points = {0};
    std::vector<float> weights = {1};
    bool weighted = true;
    bool normalized = true;
    void weight_and_normalize();
    void weight();
    void normalize();
    static float L1_dist(Hash &h1, Hash &h2);
    static float L2_dist(Hash &h1, Hash &h2);
};