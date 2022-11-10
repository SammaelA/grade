#include "model_texture_creator.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <glm/glm.hpp>
#include <vector>

ModelTex::ModelTex():
UV({"uv_coords.vs", "uv_coords.fs"}, {"in_Position", "in_Normal", "in_Tex"}), 
tex_get({"tex_from_uv.comp"},{})
{
  

  //create FBO and SSBO
  fbo = create_framebuffer();

  //bind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  //check if FBO valid
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
  
  //unbind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

ModelTex::~ModelTex()
{
  delete_framebuffer(fbo);
}
Texture ModelTex::perform_getUV(Texture &t, Texture mask, Model &m, Texture photo)
{
  getTexbyUV(t, mask, m, photo);
  return t;
}
void ModelTex::getTexbyUV(Texture &t, Texture mask, Model &m, Texture photo)
{
  //check texture type
  assert(t.type == GL_TEXTURE_2D);

  //if we don't have tmp texture or it has wrong size, we need to recreate it
  if (!tmp_tex.is_valid() || tmp_tex.type != t.type || 
      tmp_tex.get_H() != t.get_H() || tmp_tex.get_W() != t.get_W())
  {
    tmp_tex = engine::textureManager->create_texture(t.get_W(), t.get_H());
  }

  //bind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};

    Texture depthTex = engine::textureManager->create_texture(tmp_tex.get_W(), tmp_tex.get_H(), GL_DEPTH_COMPONENT16, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);
    glBindTexture(GL_TEXTURE_2D, depthTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);

  //first pass from t to tmp_tex
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex.texture, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, tmp_tex.get_W(), tmp_tex.get_H());
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glm::mat4 projection = glm::perspective(PI/2, 1.0f, 0.1f, 3000.0f);
  glm::mat4 view = glm::lookAt(glm::vec3(0, 0.5, 1.5), glm::vec3(0, 0.5, 0), glm::vec3(0, 1, 0));
  UV.use();

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(mask.type, mask.texture);
  UV.uniform("mask", 0);
  UV.uniform("projection", projection);
  UV.uniform("view", view);
  m.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  tex_get.use();
  tex_get.uniform("tex_size", glm::vec2(t.get_W(), t.get_H()));
  tex_get.texture("uv", tmp_tex);
  tex_get.texture("photo", photo);
  glBindImageTexture(2, t.texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
  
  glDispatchCompute(tmp_tex.get_W(), tmp_tex.get_H(), 1);

  glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}