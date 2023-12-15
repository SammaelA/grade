#include "render_point_cloud.h"
#include "tinyEngine/engine.h"

PointCloudRenderer::PointCloudRenderer():
shader({"point_cloud.vs", "point_cloud.fs"}, {"in_Position"})
{
  fbo = create_framebuffer();
  vao = create_vertex_array();
  glBindVertexArray(vao);
  vbo = create_buffer();
}

PointCloudRenderer::~PointCloudRenderer()
{
  delete_framebuffer(fbo);
  glDisableVertexAttribArray(vao);
  delete_vertex_array(vao);
  delete_buffer(vbo);
}

Texture PointCloudRenderer::render(const std::vector<glm::vec3> &points, const glm::mat4 &viewProj, int w, int h,
                                   glm::vec3 points_color, float points_base_size)
{
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, points.size()*sizeof(glm::vec3), points.data(), GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    Texture resTex = engine::textureManager->create_texture(w, h, GL_RGBA8);
    Texture depthTex = engine::textureManager->create_texture(w, h, GL_DEPTH_COMPONENT16, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);
    glBindTexture(GL_TEXTURE_2D, depthTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);
    glBindTexture(GL_TEXTURE_2D, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex.texture, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, resTex.texture, 0);
  glViewport(0, 0, w, h);
  glClearColor(0,0,0,1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  shader.use();
  shader.uniform("viewProj", viewProj);
  shader.uniform("points_color", points_color);
  shader.uniform("points_base_size", points_base_size);
  glDrawArrays(GL_POINTS, 0, points.size()/3);


  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

  return resTex;
}