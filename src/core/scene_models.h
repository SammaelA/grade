#pragma once
#include "common_utils/LiteMath_ext.h"
#include <vector>
#include "core/tree.h"

struct VertexData
{
    float3 pos;
    float3 normal;
    float3 tangent;
    float2 tex_coord;
};
struct SegmentVertexes
{
    float ringsize = 0;
    std::vector<VertexData> smallRing;
    std::vector<VertexData> bigRing;
    Segment s;
};