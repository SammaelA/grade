#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include <vector>

class ModelTex
{
public:
  ModelTex();
  ~ModelTex();
  Texture getTexbyUV(Texture mask, Model &m, Texture photo, int overdraw, const CameraSettings &camera);
  typedef struct tex_data
  {
    double w0; 
    double h0; 
    double w1; 
    double h1; 
    int x_sym; 
    int y_sym;
  } tex_data;
  Texture symTexComplement(Texture tex, std::vector<tex_data> texs_data);
private:
  Texture tmp_tex;
  GLuint fbo;
  Shader UV;
  Shader tex_get;
  PostFx cpy1;
  PostFx cpy2;
  PostFx texture_postprocess;
  PostFx texture_mirror;
  PostFx tex_com;
  PostFx texs_div;
  PostFx tex_to_mask;
  PostFx photo_transform;
};