#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <string>

namespace dgen
{
  std::vector<float> load_obj(const std::string &filename);
  void save_obj(const std::string &filename, const std::vector<float> &model);
};