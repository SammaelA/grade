#pragma once
#include "common_utils/bbox.h"
#include "core/scene.h"

namespace SceneGenHelper
{
  uint64_t pack_id(unsigned _empty, unsigned category, unsigned type, unsigned id);
  void unpack_id(uint64_t packed_id, unsigned &_empty, unsigned &category, unsigned &type, unsigned &id);
  void get_AABB_list_from_instance(ComplexModel &m, glm::mat4 transform, std::vector<AABB> &boxes, float column_size, 
                                   float inflation_q = 1.0f);

  bool is_terrain(glm::vec4 world_pos_type);
};