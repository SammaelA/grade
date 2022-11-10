#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

class ModelTex
{
public:
  ModelTex();
  ~ModelTex();
  void perform_getUV_inplace(Texture &t, Texture mask, Model &m, Texture photo);
  Texture perform_getUV(Texture &t, Texture mask, Model &m, Texture photo);
private:
  void getTexbyUV(Texture &t, Texture mask, Model &m, Texture photo);
  Texture tmp_tex;
  GLuint fbo;
  Shader UV;
  Shader tex_get;
};