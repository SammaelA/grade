#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

struct ImgExp
{
  static Texture ImgExpanding(Texture image, int res_size);//returns prepocessed silhouette and reference texture
};