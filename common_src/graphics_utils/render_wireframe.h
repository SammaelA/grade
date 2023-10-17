#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include <vector>

class WireframeRenderer
{
public:
  WireframeRenderer();
  ~WireframeRenderer();
  Texture render(Model &m, const glm::mat4 &viewProj, int tex_w, int tex_h);
private:
  GLuint fbo;
  Shader wireframe_shader;
};