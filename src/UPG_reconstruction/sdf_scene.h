#pragma once
#include "upg.h"
#include "common_utils/bbox.h"

namespace upg
{
  struct SdfObject
  {
    unsigned type;
    unsigned params_offset; //in parameters vector
    unsigned params_count;
    float distance_mult = 1.0f;
    float distance_add = 0.0f;
    AABB bbox;
    glm::mat4 transform;
    bool complement = false;
  };
  struct SdfConjunction
  {
    unsigned offset; //in objects vector
    unsigned size;
  };
  struct SdfScene
  {
    std::vector<float> parameters;
    std::vector<SdfObject> objects;
    std::vector<SdfConjunction> conjunctions;
  };

  SdfScene create_sdf_scene(const UPGStructure &structure, const UPGParametersRaw &params);
}