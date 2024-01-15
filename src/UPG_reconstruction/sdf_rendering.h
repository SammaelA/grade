#pragma once
#include "upg.h"
#include "sdf_node.h"
#include "tinyEngine/camera.h"
#include "common_utils/bbox.h"

namespace upg
{
  AABB get_point_cloud_bbox(const std::vector<glm::vec3> &points);
  void sdf_to_point_cloud(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                          std::vector<glm::vec3> *outside_points = nullptr);
  void sdf_to_point_cloud_with_dist(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                                    std::vector<float> *distances);
  Texture render_sdf(const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, bool lambert = true);
  float get_sdf_image_based_quality(ProceduralSdf reference_sdf, ProceduralSdf sdf);
}