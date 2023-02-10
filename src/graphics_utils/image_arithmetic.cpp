#include "image_arithmetic.h"
#include "tinyEngine/engine.h"

#define ADD 0
#define MUL 1
#define DIV 2

void image_op(Texture &I, Texture &I_1, Texture &I_2, float a, float b, int op)
{
  assert(I_1.type == GL_TEXTURE_2D);
  assert(I_2.type == GL_TEXTURE_2D);
  assert(I.type == GL_TEXTURE_2D);
  PostFx filter("image_arithmetics.fs");
  float w = I.get_W();
  float h = I.get_H();

  GLuint fbo = create_framebuffer();
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, I.texture, 0);
  glViewport(0, 0, w, h);
  filter.use();
  filter.get_shader().texture("tex_1", I_1);
  filter.get_shader().texture("tex_2", I_2);
  filter.get_shader().uniform("op", op);
  filter.get_shader().uniform("a", a);
  filter.get_shader().uniform("b", b);
  filter.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  delete_framebuffer(fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

void ImageArithmetics::add(Texture &I, Texture &I_1, Texture &I_2, float a, float b)
{
  return image_op(I, I_1, I_2, a, b, ADD);
}

void ImageArithmetics::mul(Texture &I, Texture &I_1, Texture &I_2, float a)
{
  return image_op(I, I_1, I_2, a, 1, MUL);
}

void ImageArithmetics::div(Texture &I, Texture &I_1, Texture &I_2, float a)
{
  return image_op(I, I_1, I_2, a, 1, DIV);
}