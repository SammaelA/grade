#pragma once

#include "common_utils/bbox.h"
#include <vector>
#include <functional>

namespace upg
{
std::vector<AABB> get_bbox_list(std::function<float(const float3 &)> sdf,
                                const AABB &sdf_bbox, int bbox_count, int cuts_quant);
}