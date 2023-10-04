#include "bilateral_filter.h"
#include "tinyEngine/engine.h"

Texture BilateralFilter::perform(Texture &t, float sigma_d, float sigma_r)
{
  assert(t.type == GL_TEXTURE_2D);
  PostFx filter("bilateral_filter.fs");
  float w = t.get_W();
  float h = t.get_H();

  Texture res_tex = engine::textureManager->create_texture(w, h);
  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res_tex.texture, 0);
  glViewport(0, 0, w, h);
  filter.use();
  filter.get_shader().texture("tex",t);
  filter.get_shader().uniform("tex_size", glm::vec2(w, h));
  filter.get_shader().uniform("sigma_d",sigma_d);
  filter.get_shader().uniform("sigma_r",sigma_r);
  filter.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  delete_framebuffer(fbo);
  return res_tex;
}