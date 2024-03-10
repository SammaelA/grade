#pragma once
#include <vector>
#include "common_utils/LiteMath_ext.h"
#include <string>

#define FLOAT_PER_VERTEX (3+3+2)

namespace dgen
{
  std::vector<float> load_obj(const std::string &filename);
  void set_face_normals(std::vector<float> &model);
  void save_obj(const std::string &filename, const std::vector<float> &model);
};