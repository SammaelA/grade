#pragma once
#include <string>
#include "common_utils/LiteMath_ext.h"

namespace compare_utils
{
  float loss_psnr(const std::string &img1, const std::string &img2);
  float loss_silhouette_psnr(const std::string &img1, const std::string &img2, float3 background_color);
  float loss_silhouette_iou(const std::string &img1, const std::string &img2, float3 background_color);

  void turntable_loss(const std::string &path1, const std::string &path2, int image_count);
} // namespace compare_utils
