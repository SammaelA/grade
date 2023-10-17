#include "render_wireframe.h"
#include "tinyEngine/engine.h"

WireframeRenderer::WireframeRenderer():
wireframe_shader({"wireframe.vs", "wireframe.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
  fbo = create_framebuffer();
}

WireframeRenderer::~WireframeRenderer()
{
  delete_framebuffer(fbo);
}

Texture WireframeRenderer::render(Model &m, const glm::mat4 &viewProj, int w, int h)
{
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
  glClearColor(1,1,1,1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  wireframe_shader.use();
  wireframe_shader.uniform("viewProj", viewProj);

  wireframe_shader.uniform("wireframe_color", glm::vec3(1,1,1));
  m.render();
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  wireframe_shader.uniform("wireframe_color", glm::vec3(0,0,0));
  m.render();

  glClearColor(0,0,0,0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

  return resTex;
}