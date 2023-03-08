#include "image_expand.h"
#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <glm/glm.hpp>
#include <vector>

Texture ImgExp::ImgExpanding(Texture image, int res_size, float color_thr, float blur_sigma)
{
  GLuint fbo;
  Shader get_borders({"take_mask_borders.comp"},{});
  PostFx img_exp("restrict_tex.fs");

  //create FBO and SSBO
  fbo = create_framebuffer();

  //bind FBO
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  //check if FBO valid
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));

  float w = image.get_W();
  float h = image.get_H();

  SilhouetteExtractor se = SilhouetteExtractor(blur_sigma >=0 ? blur_sigma : MAX(1, MIN(w, h)/256.0), 0, 0, color_thr);
  Texture mask = se.get_silhouette_simple(image, w, h);
  engine::textureManager->save_png(mask, "ie_mask");
  engine::textureManager->save_png(image, "ie_image");

    float borderColorDepth[] = {0.0f, 0.0f, 0.0f, 0.0f};
    glBindTexture(GL_TEXTURE_2D, mask.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);
    glBindTexture(GL_TEXTURE_2D, 0);

  uint ssbo = 0;
  glGenBuffers(1, &ssbo);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
  get_borders.use();
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
  glBufferData(GL_SHADER_STORAGE_BUFFER, 16, NULL, GL_DYNAMIC_READ);
  get_borders.uniform("size", glm::vec2(w, h));
  get_borders.texture("mask", mask);
  
  glDispatchCompute(1, 1, 1);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  int arr[4] = {-2, -2, -2, -2};
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
  glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, 16, arr);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
  glDeleteBuffers(1, &ssbo);

  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  PostFx copy("copy.fs");
  Texture tmp_tex = engine::textureManager->create_texture(res_size, res_size);

  float real_w = (arr[1] - arr[0]);
  float real_h = (arr[3] - arr[2]);
  glm::vec4 tex_transform(arr[0]/(float)w, arr[2]/(float)h, real_w/w, real_h/h);

  float border = 0.05;
  float sz_x = real_w > real_h ? 1 : real_w/real_h;
  float sz_y = real_w > real_h ? real_h/real_w : 1;

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, tmp_tex.get_W(), tmp_tex.get_H());
  glClearColor(0, 0, 0, 1);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  glViewport((0.5*(1 - sz_x) + border)*tmp_tex.get_W(), (0.5*(1-sz_y) + border)*tmp_tex.get_H(), 
             (sz_x - 2*border)*tmp_tex.get_W(), (sz_y - 2*border)*tmp_tex.get_H());

  img_exp.use();
  img_exp.get_shader().texture("tex", image);
  img_exp.get_shader().texture("mask", mask);
  img_exp.get_shader().uniform("tex_transform", tex_transform);
  img_exp.render();

  delete_framebuffer(fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
  engine::textureManager->save_png(tmp_tex, "ie_rsult");
  return tmp_tex;
}