#include "model_texture_creator.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
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
Texture ModelTex::perform_getUV(Texture &t, Texture mask, Model m, Texture photo)
{
  getTexbyUV(t, mask, m, photo);
  return t;
}
void ModelTex::getTexbyUV(Texture &t, Texture mask, Model m, Texture photo)
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

  //first pass from t to tmp_tex
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  
  UV.use();

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(mask.type, mask.texture);
  UV.uniform("mask", 0);
  
  m.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  tex_get.use();

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(tmp_tex.type, tmp_tex.texture);
  UV.uniform("uv", 1);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(photo.type, photo.texture);
  UV.uniform("photo", 1);

  glActiveTexture(GL_TEXTURE2);
  glBindTexture(t.type, t.texture);
  UV.uniform("tex", 2);
  
  glDispatchCompute(tmp_tex.get_W(), tmp_tex.get_H(), 1);

  glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}