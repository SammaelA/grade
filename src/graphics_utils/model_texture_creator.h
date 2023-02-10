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
  Texture getTexbyUV(Texture mask, Model &m, Texture photo, const CameraSettings &camera, Texture &res_mask,
                     int res_overdraw = 1, int reconstruction_overdraw = 3);
  typedef struct tex_data
  {
    double w0; 
    double h0; 
    double w1; 
    double h1; 
    int x_sym; 
    int y_sym;
  } tex_data;
  Texture symTexComplement(Texture tex, Texture mask, std::vector<tex_data> texs_data);
private:
  GLuint fbo;
  Shader UV;
  Shader tex_get;
  PostFx cpy1;
  PostFx cpy2;
  PostFx texture_postprocess;
  PostFx tex_com;
  PostFx texs_div;
  PostFx photo_transform;
};