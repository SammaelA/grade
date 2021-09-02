#pragma once
#include <vector>

struct Hash
{
    std::vector<float> data;
    std::vector<int> start_points = {0};
    std::vector<float> weights = {1};
    bool weighted = false;
    bool normalized = false;
    void weight_and_normalize();
    static float L1_dist(Hash &h1, Hash &h2);
    static float L2_dist(Hash &h1, Hash &h2);
};