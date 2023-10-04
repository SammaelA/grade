#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

struct ImgExp
{
  static Texture ImgExpanding(Texture image, int res_size, float color_thr = 0.25, float blur_sigma = -1);//returns prepocessed silhouette and reference texture
};