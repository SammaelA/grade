#include "image_expand.h"
#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <glm/glm.hpp>
#include <vector>

ImgExp::ImgExp():
get_borders({"take_mask_borders.comp"},{}),
img_exp("restrict_tex.fs")
{
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
  
  //unbind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

ImgExp::~ImgExp()
{
  delete_framebuffer(fbo);
}

Texture ImgExp::ImgExpanding(Texture image)
{
  //check texture type
  logerr("%s", glGetString(GL_VERSION));
  Texture tmp = engine::textureManager->create_texture(image.get_W(), image.get_H());

  float w = image.get_W();
  float h = image.get_H();
  //bind FBO
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
  Texture mask = se.get_silhouette(image, MAX(w, h), MAX(w, h));

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
  get_borders.uniform("size", glm::vec2(MAX(w, h), MAX(w, h)));
  get_borders.texture("mask", mask);
  
  glDispatchCompute(1, 1, 1);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  int arr[4] = {-2, -2, -2, -2};
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
  glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, 16, arr);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

  float W = (arr[1] - arr[0]) / w;
  float H = (arr[3] - arr[2]) / h;
  float P = MAX(W, H);
  //std::cout << arr[0] << " " << arr[1] << " " << arr[2] << " " << arr[3] << " | " << W << " " << H << " " << P << " | ";
  W = ((arr[1] + arr[0] + w - MAX(w, h)) / w - P) / 2;
  H = ((arr[2] + arr[3] + h - MAX(w, h)) / h - P) / 2;
  //std::cout << W << " " << H;
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glBindTexture(GL_TEXTURE_2D, image.texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glBindTexture(GL_TEXTURE_2D, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp.texture, 0);
    glViewport(0, 0, w, h);
    img_exp.use();
    img_exp.get_shader().uniform("x_sh", W);
    img_exp.get_shader().uniform("y_sh", H);
    img_exp.get_shader().uniform("x_sz", P);
    img_exp.get_shader().uniform("y_sz", P);
    img_exp.get_shader().texture("tex", image);
    img_exp.render();

  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

  return tmp;
}