#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

class GaussFilter
{
public:
  GaussFilter(float sigma, float precision = 0.003);
  ~GaussFilter();
  void perform_gauss_blur(Texture &t);
private:
  std::vector<float> kernel;//1D gauss kernel
  Texture tmp_tex;//intermediate texture
  GLuint fbo;//frame buffer, tmp_tex and buffer will be attached to it
  GLuint kernel_buf;//buffer where kernel values will be stored
  PostFx blur;//postfx that contains blur shader (it applies shader to a full-screen quad with texture)
};