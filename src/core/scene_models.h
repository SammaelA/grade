#pragma once
#include "common_utils/LiteMath_ext.h"
#include <vector>
#include "core/tree.h"

struct VertexData
{
    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec3 tangent;
    glm::vec2 tex_coord;
};
struct SegmentVertexes
{
    float ringsize = 0;
    std::vector<VertexData> smallRing;
    std::vector<VertexData> bigRing;
    Segment s;
};