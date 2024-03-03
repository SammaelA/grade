#pragma once
#include "LiteMath_ext.h"
#include <string>
#include <vector>

  struct SdfObject
  {
    unsigned type;
    unsigned params_offset; //in parameters vector
    unsigned params_count;
    float distance_mult = 1.0f;
    float distance_add = 0.0f;
    LiteMath::AABB bbox;
    LiteMath::float4x4 transform;
    bool complement = false;
  };
  struct SdfConjunction
  {
    unsigned offset; //in objects vector
    unsigned size;
    LiteMath::AABB bbox;
  };
  struct SdfScene
  {
    std::vector<float> parameters;
    std::vector<SdfObject> objects;
    std::vector<SdfConjunction> conjunctions;
  };

  //evaluating and rendering
  float dist_prim(const SdfScene &sdf, const SdfObject &prim, LiteMath::float3 p);
  float get_dist(const SdfScene &sdf, LiteMath::float3 p);
  bool sdf_sphere_tracing(const SdfScene &sdf, const LiteMath::AABB &sdf_bbox, const LiteMath::float3 &pos, const LiteMath::float3 &dir, 
                          LiteMath::float3 *surface_pos = nullptr);

  //save/load scene
  void save_sdf_scene(const SdfScene &scene, const std::string &path);
  void load_sdf_scene(SdfScene &scene, const std::string &path);