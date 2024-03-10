#pragma once
#include <vector>
#include "common_utils/LiteMath_ext.h"
#include "common_utils/bbox.h"
#include "common_utils/blk.h"
#include "tinyEngine/camera.h"

namespace visualizer
{
  void shift(std::vector<float> &model, glm::vec3 shift);
  void scale(std::vector<float> &model, glm::vec3 scale);
  void transform(std::vector<float> &model, glm::mat4 transform_mat);
  
  AABB get_bbox(const std::vector<float> &model);
  glm::vec3 get_center_of_mass_vertices(const std::vector<float> &model);
  void normalize_model(std::vector<float> &model);

  void save_camera_settings(const CameraSettings &camera, Block &blk);
  CameraSettings load_camera_settings(const Block &blk);
};