#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include <vector>

SilhouetteExtractor::SilhouetteExtractor(float _blur_sigma, float low_thr, float high_thr, float _color_thr):
blur_sigma(_blur_sigma),
color_thr(_color_thr),
gauss(MAX(0.1, _blur_sigma)),
detect_object("detect_object.fs"),
blur_mask_edges("blur_mask_edges.fs"),
remove_holes("remove_holes.fs")
{

}

Texture SilhouetteExtractor::get_silhouette(Texture &t, int res_w, int res_h)
{
  assert(t.type == GL_TEXTURE_2D);
  float w = t.get_W();
  float h = t.get_H();

  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  Texture tmp_tex = engine::textureManager->create_texture(res_w, res_h);
  Texture res = engine::textureManager->create_texture(res_w, res_h);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, res_w, res_h);
  glClearColor(0, 0, 0, 1);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  detect_object.use();
  detect_object.get_shader().texture("tex_color", t);
  detect_object.get_shader().uniform("tex_size",float2(w, h));
  detect_object.get_shader().uniform("color_thr", color_thr);
  detect_object.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
  remove_holes.use();
  remove_holes.get_shader().texture("tex_mask", tmp_tex);
  remove_holes.get_shader().uniform("tex_size",float2(res_w, res_h));
  remove_holes.get_shader().uniform("search_radius", (int)MIN(0.03*res_w, 0.03*res_h));
  remove_holes.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  remove_holes.use();
  remove_holes.get_shader().texture("tex_mask", res);
  remove_holes.get_shader().uniform("tex_size",float2(res_w, res_h));
  remove_holes.get_shader().uniform("search_radius", (int)MIN(0.03*res_w, 0.03*res_h));
  remove_holes.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
  blur_mask_edges.use();
  blur_mask_edges.get_shader().texture("tex_mask", tmp_tex);
  blur_mask_edges.get_shader().uniform("tex_size",float2(res_w, res_h));
  blur_mask_edges.get_shader().uniform("search_radius", -2);
  blur_mask_edges.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  if (blur_sigma > 0.1)
    gauss.perform_gauss_blur_inplace(res);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  delete_framebuffer(fbo);
  engine::view->next_frame();

  return res;
}