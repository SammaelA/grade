#pragma once

struct Texture;
class ImageMetric
{
public:
  enum Metric
  {
    MAE,
    MSE,
    PSNR,
    IOU
  };

  static float get(const Texture &t1, const Texture &t2, Metric metric = MSE);
};