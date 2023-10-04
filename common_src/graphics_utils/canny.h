#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "gauss_blur_precise.h"
#include <vector>

class Canny
{
public:
  Canny(float blur_sigma = 1.0, float low_thr = 0.125, float high_thr = 0.375);
  Texture detect_edges(Texture &input);
private:
  float low_thr, high_thr;
  GaussFilter gauss;
  PostFx get_grad_discrete, non_maximum_supression, hysteresis;
};