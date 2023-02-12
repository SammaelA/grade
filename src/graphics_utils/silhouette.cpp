#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <vector>

SilhouetteExtractor::SilhouetteExtractor(float blur_sigma, float low_thr, float high_thr):
gauss(7.5f),
canny(blur_sigma, low_thr, high_thr),
fill_edges("fill_edges.fs"),
fill_silhouette("fill_edges.fs"),
detect_object("detect_object.fs"),
remove_holes("remove_holes.fs"),
copy("copy.fs"), 
metric("ref_image_preprocess_metric.fs")
{

}
Texture SilhouetteExtractor::get_silhouette(Texture &t, int res_w, int res_h)
{
  engine::textureManager->save_png(t, "initial_tex0");
  assert(t.type == GL_TEXTURE_2D);
  float w = t.get_W();
  float h = t.get_H();

  Texture tmp_tex = engine::textureManager->create_texture(w, h);
  Texture canny_tex = canny.detect_edges(t);
  
  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, w, h);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  fill_edges.use();
  fill_edges.get_shader().texture("tex", canny_tex);
  fill_edges.get_shader().uniform("tex_size",glm::vec2(w, h));
  fill_edges.get_shader().uniform("search_radius",4);
  fill_edges.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, canny_tex.texture, 0);
  detect_object.use();
  detect_object.get_shader().texture("tex_edges", tmp_tex);
  detect_object.get_shader().texture("tex_color", t);
  detect_object.get_shader().uniform("tex_size",glm::vec2(w, h));
  detect_object.get_shader().uniform("color_thr",0.03f);
  detect_object.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  remove_holes.use();
  remove_holes.get_shader().texture("tex_mask", canny_tex);
  remove_holes.get_shader().uniform("tex_size",glm::vec2(w, h));
  remove_holes.get_shader().uniform("search_radius", 3);
  remove_holes.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  Texture res = engine::textureManager->create_texture(res_w, res_h);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
  glViewport(0, 0, res_w, res_h);
  glClearColor(0, 0, 0, 1);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  {
    float mul_x = MIN(1, w/h);
    float mul_y = MIN(1, h/w);
    glViewport(0.5*res_w*(1-mul_x), 0.5*res_w*(1-mul_y), mul_x*res_w, mul_y*res_h);
  }
  copy.use();
  copy.get_shader().texture("tex",tmp_tex);
  copy.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  engine::textureManager->save_png(tmp_tex, "fill_edges");
  engine::textureManager->save_png(t, "initial_tex");
  engine::textureManager->save_png(canny_tex, "initial_mask");
  engine::textureManager->save_png(res, "res_tex");

  Texture result = engine::textureManager->create_texture(res_w, res_h);

  glBindTexture(GL_TEXTURE_2D, res.texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  float clr[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, clr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
  glBindTexture(GL_TEXTURE_2D, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, result.texture, 0);
  glViewport(0, 0, res_w, res_h);
  metric.use();
  metric.get_shader().uniform("max_rad", 32.0);
  metric.get_shader().uniform("size_x", (float)res_w);
  metric.get_shader().uniform("size_y", (float)res_h);
  metric.get_shader().texture("sil", res);
  metric.render();
  
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  delete_framebuffer(fbo);
  engine::view->next_frame();

  return result;
}