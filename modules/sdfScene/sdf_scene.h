#pragma once
#include "common_utils/bbox.h"
#include "LiteMath/LiteMath.h"
#include <string>
#include <vector>

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
    AABB bbox;
  };
  struct SdfScene
  {
    std::vector<float> parameters;
    std::vector<SdfObject> objects;
    std::vector<SdfConjunction> conjunctions;
  };

  //evaluating and rendering
  float dist_prim(const SdfScene &sdf, const SdfObject &prim, glm::vec3 p);
  float get_dist(const SdfScene &sdf, glm::vec3 p);
  bool sdf_sphere_tracing(const SdfScene &sdf, const AABB &sdf_bbox, const glm::vec3 &pos, const glm::vec3 &dir, 
                          glm::vec3 *surface_pos = nullptr);

  //save/load scene
  void save_sdf_scene(const SdfScene &scene, const std::string &path);
  void load_sdf_scene(SdfScene &scene, const std::string &path);