#pragma once
#include "sdfScene/LiteMath_ext.h"
#include "common_utils/bbox.h"
#include "common_utils/LiteMath_ext.h"

static glm::vec3 conv(const LiteMath::float3 &v)
{
  return glm::vec3(v.x, v.y, v.z);
}
static LiteMath::float3 conv(const glm::vec3 &v)
{
  return LiteMath::float3(v.x, v.y, v.z);
}
static AABB conv(const LiteMath::AABB &box)
{
  return AABB(conv(box.min_pos), conv(box.max_pos));
}
static LiteMath::AABB conv(const AABB &box)
{
  return LiteMath::AABB(conv(box.min_pos), conv(box.max_pos));
}