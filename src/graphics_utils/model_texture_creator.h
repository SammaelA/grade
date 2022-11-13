#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include <vector>

class ModelTex
{
public:
  ModelTex();
  ~ModelTex();
  Texture getTexbyUV(Texture mask, Model &m, Texture photo, int overdraw);
private:
  Texture tmp_tex;
  GLuint fbo;
  Shader UV;
  Shader tex_get;
  PostFx texture_postprocess;
  PostFx texture_mirror;
  PostFx photo_transform;
};