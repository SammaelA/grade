#include "resize_image.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/postfx.h"

Texture ImageResizer::resize(Texture t, int new_w, int new_h, Type type, float4 base_color)
{
  assert(t.type == GL_TEXTURE_2D);
  float w = t.get_W();
  float h = t.get_H();

  Texture tmp_tex = engine::textureManager->create_texture(new_w, new_h);
  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, new_w, new_h);
  glClearColor(base_color.x, base_color.y, base_color.z, base_color.w);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  if (type == Type::STRETCH)
  {
    glViewport(0, 0, new_w, new_h);
  }
  else if (type == Type::CENTERED)
  {
    float mul_x = MIN(1, w/h);
    float mul_y = MIN(1, h/w);
    glViewport(0.5*new_w*(1-mul_x), 0.5*new_w*(1-mul_y), mul_x*new_w, mul_y*new_h);
  }
  PostFx copy("copy.fs");
  copy.use();
  copy.get_shader().texture("tex",t);
  copy.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  delete_framebuffer(fbo);
  return tmp_tex;
}