#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include <vector>

class DepthLossCalculator
{
public:
  DepthLossCalculator();
  ~DepthLossCalculator();
  Texture get_depth(Model &m, const CameraSettings &camera, int tex_size_x, int tex_size_y);
  float get_loss(Model &m, Texture reference_depth, const CameraSettings &camera);
private:
  GLuint fbo;
  Shader render_model;
  PostFx depth_postprocess;
};