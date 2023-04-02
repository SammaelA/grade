#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "canny.h"
#include <vector>

class SilhouetteExtractor
{
public:
  SilhouetteExtractor(float blur_sigma = 1.0, float low_thr = 0.125, float high_thr = 0.375, float color_thr = 0.05);
  Texture get_silhouette(Texture &input, int w, int h);
private:
  float blur_sigma, color_thr;
  GaussFilter gauss;
  PostFx detect_object, blur_mask_edges, remove_holes;
};