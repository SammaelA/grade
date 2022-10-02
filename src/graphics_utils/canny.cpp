#include "canny.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <vector>

Canny::Canny(float blur_sigma, float _low_thr, float _high_thr):
low_thr(_low_thr),
high_thr(_high_thr),
gauss(blur_sigma),
get_grad_discrete("get_grad_discrete.fs"),
non_maximum_supression("non_maximum_supression.fs"),
hysteresis("hysteresis.fs")
{

}

Texture Canny::detect_edges(Texture &t)
{
  assert(t.type == GL_TEXTURE_2D);
  float w = t.get_W();
  float h = t.get_H();

  Texture grad_tex = engine::textureManager->create_texture(w, h);
  Texture res_tex = gauss.perform_gauss_blur(t);
  
  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, grad_tex.texture, 0);
  glViewport(0, 0, w, h);

  get_grad_discrete.use();
  get_grad_discrete.get_shader().texture("tex", res_tex);
  get_grad_discrete.get_shader().uniform("tex_size",glm::vec2(w, h));
  get_grad_discrete.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res_tex.texture, 0);
  non_maximum_supression.use();
  non_maximum_supression.get_shader().texture("tex", grad_tex);
  non_maximum_supression.get_shader().uniform("tex_size",glm::vec2(w, h));
  non_maximum_supression.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, grad_tex.texture, 0);
  hysteresis.use();
  hysteresis.get_shader().texture("tex", res_tex);
  hysteresis.get_shader().uniform("tex_size",glm::vec2(w, h));
  hysteresis.get_shader().uniform("thresholds",glm::vec2(low_thr, high_thr));
  hysteresis.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  delete_framebuffer(fbo);

  return grad_tex;
}