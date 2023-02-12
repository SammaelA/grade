#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

class ImgExp
{
public:
  ImgExp();
  ~ImgExp();
  Texture ImgExpanding(Texture image);
private:
  GLuint fbo;
  Shader get_borders;
  PostFx img_exp;
};