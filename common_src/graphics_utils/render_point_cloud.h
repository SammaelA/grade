#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include <vector>

class PointCloudRenderer
{
public:
  PointCloudRenderer();
  ~PointCloudRenderer();
  Texture render(const std::vector<float3> &points, const float4x4 &viewProj, int tex_w, int tex_h,
                 float3 points_color = float3(1,1,1), float points_base_size = 0.1);
private:
  GLuint fbo, vao, vbo;
  Shader shader;
};