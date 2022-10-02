#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <vector>

SilhouetteExtractor::SilhouetteExtractor(float blur_sigma, float low_thr, float high_thr):
canny(blur_sigma, low_thr, high_thr),
fill_edges("fill_edges.fs"),
fill_silhouette("fill_edges.fs"),
detect_object("detect_object.fs")
{

}
Texture SilhouetteExtractor::get_silhouette(Texture &t)
{
  engine::textureManager->save_png(t, "initial_tex0");
  assert(t.type == GL_TEXTURE_2D);
  float w = t.get_W();
  float h = t.get_H();

  Texture tmp_tex = engine::textureManager->create_texture(w, h);
  Texture res = canny.detect_edges(t);
  
  GLuint fbo = create_framebuffer();
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, w, h);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  fill_edges.use();
  fill_edges.get_shader().texture("tex", res);
  fill_edges.get_shader().uniform("tex_size",glm::vec2(w, h));
  fill_edges.get_shader().uniform("search_radius",4);
  fill_edges.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
  detect_object.use();
  detect_object.get_shader().texture("tex_edges", tmp_tex);
  detect_object.get_shader().texture("tex_color", t);
  detect_object.get_shader().uniform("tex_size",glm::vec2(w, h));
  detect_object.get_shader().uniform("color_thr",0.01f);
  detect_object.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  engine::textureManager->save_png(tmp_tex, "fill_edges");
  engine::textureManager->save_png(t, "initial_tex");
  engine::textureManager->save_png(res, "initial_mask");

  delete_framebuffer(fbo);

  return res;
}