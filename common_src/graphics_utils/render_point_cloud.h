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
  Texture render(const std::vector<glm::vec3> &points, const glm::mat4 &viewProj, int tex_w, int tex_h,
                 glm::vec3 points_color = glm::vec3(1,1,1), float points_base_size = 0.1);
private:
  GLuint fbo, vao, vbo;
  Shader shader;
};