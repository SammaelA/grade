#include "model_texture_creator.h"
#include "gauss_blur_precise.h"
#include "silhouette.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include "resize_image.h"
#include <glm/glm.hpp>
#include <vector>

ModelTex::ModelTex():
UV({"uv_coords.vs", "uv_coords.fs"}, {"in_Position", "in_Normal", "in_Tex"}), 
tex_get({"tex_from_uv.comp"},{}),
photo_transform("copy.fs"),
texture_postprocess("texture_postprocess.fs"),
tex_com("tex_com_2.fs"),
texs_div("tex_div_tex.fs"),
copy("copy.fs")
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

ModelTex::~ModelTex()
{
  delete_framebuffer(fbo);
}

Texture ModelTex::getTexbyUV(Texture mask, Model &m, Texture photo, const CameraSettings &camera, Texture &res_mask,
                             int res_od, int rec_od)
{
  //check texture type
  Texture t = engine::textureManager->create_texture(res_od*rec_od*photo.get_W(), res_od*rec_od*photo.get_H());
  Texture tmp_mask = engine::textureManager->create_texture(res_od*rec_od*photo.get_W(), res_od*rec_od*photo.get_H());

  int w = t.get_W();
  int h = t.get_H();
  //bind FBO
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    Texture photo_transformed = engine::textureManager->create_texture(w, h);
    Texture UV_tex = engine::textureManager->create_texture(w, h, GL_RGB32F);
    Texture depthTex = engine::textureManager->create_texture(w, h, GL_DEPTH_COMPONENT16, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);
    glBindTexture(GL_TEXTURE_2D, depthTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);
    glBindTexture(GL_TEXTURE_2D, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex.texture, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, UV_tex.texture, 0);
  glViewport(0, 0, w, h);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glm::mat4 projection = glm::perspective(camera.fov_rad, 1.0f, camera.z_near, camera.z_far);
  glm::mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);

  UV.use();
  UV.uniform("projection", projection);
  UV.uniform("view", view);
  m.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, photo_transformed.texture, 0);
  glDisable(GL_DEPTH_TEST);
  photo_transform.use();
  photo_transform.get_shader().uniform("tex_transform", glm::vec4(0, 1, 1, -1));
  photo_transform.get_shader().texture("tex", photo);
  photo_transform.render();
  
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  tex_get.use();
  tex_get.uniform("tex_size", glm::vec2(w, h));
  tex_get.texture("uv", UV_tex);
  tex_get.texture("photo", photo_transformed);
  glBindImageTexture(2, t.texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
  glBindImageTexture(3, tmp_mask.texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
  
  glDispatchCompute(w, h, 1);

  glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

  Texture tmp_tex1 = engine::textureManager->create_texture(w, h);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex1.texture, 0);
  glViewport(0, 0, tmp_tex1.get_W(), tmp_tex1.get_H());
  texture_postprocess.use();
  texture_postprocess.get_shader().texture("tex", t);
  texture_postprocess.get_shader().uniform("tex_size", glm::vec2(w, h));
  texture_postprocess.get_shader().uniform("radius", rec_od + 1);
  texture_postprocess.get_shader().uniform("alpha_thr", 0.1);
  texture_postprocess.render();
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  Texture tmp_tex2 = engine::textureManager->create_texture(res_od*photo.get_W(), res_od*photo.get_H());
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex2.texture, 0);
  glViewport(0, 0, tmp_tex2.get_W(), tmp_tex2.get_H());
  texture_postprocess.use();
  texture_postprocess.get_shader().texture("tex", tmp_tex1);
  texture_postprocess.get_shader().uniform("tex_size", glm::vec2(w, h));
  texture_postprocess.get_shader().uniform("radius", rec_od + 1);
  texture_postprocess.get_shader().uniform("alpha_thr", 0.3);
  texture_postprocess.render();

  res_mask = engine::textureManager->create_texture(res_od*photo.get_W(), res_od*photo.get_H());
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res_mask.texture, 0);
  glViewport(0, 0, res_mask.get_W(), res_mask.get_H());
  PostFx average("average_filter.fs");
  average.use();
  average.get_shader().texture("tex", tmp_mask);
  average.get_shader().uniform("radius", rec_od + 1);
  average.get_shader().uniform("tex_size", glm::vec2(tmp_mask.get_W(), tmp_mask.get_H()));
  average.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glEnable(GL_DEPTH_TEST);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

  return tmp_tex2;
}

Texture ModelTex::symTexComplement(Texture tex, Texture mask, std::vector<tex_data> texs_data)
{
  Texture t = engine::textureManager->create_texture(tex.get_W(), tex.get_H());
  int W = tex.get_W();
  int H = tex.get_H();

  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glDisable(GL_DEPTH_TEST);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.texture, 0);
  glViewport(0, 0, W, H);
  copy.use();
  copy.get_shader().texture("tex", tex);
  copy.render();
  
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  for (auto it : texs_data)
  {
    double w = it.w1 - it.w0;
    double h = it.h1 - it.h0;
    if (w < 0 || h < 0) continue;
    Texture tmp_tex = engine::textureManager->create_texture(w * W, h * H);
    Texture sm_tex = engine::textureManager->create_texture(w * W, h * H);

    glBindTexture(GL_TEXTURE_2D, t.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (it.x_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    }
    if (it.y_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sm_tex.texture, 0);
    glViewport(0, 0, W * w, H * h);
    copy.use();
    copy.get_shader().uniform("tex_transform", glm::vec4(it.w0, it.h0, w, h));
    copy.get_shader().texture("tex", t);
    copy.render();
    
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
//////

    Texture t1 = engine::textureManager->create_texture(W * w, H * h);

    Texture mask1 = engine::textureManager->create_texture(W * w, H * h);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, mask1.texture, 0);
    glViewport(0, 0, W * w, H * h);
    copy.use();
    copy.get_shader().uniform("tex_transform", glm::vec4(it.w0, it.h0, w, h));
    copy.get_shader().texture("tex", mask);
    copy.render();

    GaussFilter gaussian(3);

    Texture t2 = gaussian.perform_gauss_blur(sm_tex);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t1.texture, 0);
    texs_div.use();
    texs_div.get_shader().texture("tex1", sm_tex);
    texs_div.get_shader().texture("tex2", t2);
    texs_div.render();
    
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    Texture mask2 = gaussian.perform_gauss_blur(mask1);

    glBindTexture(GL_TEXTURE_2D, t1.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (it.x_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    }
    if (it.y_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    glBindTexture(GL_TEXTURE_2D, t2.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (it.x_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    }
    if (it.y_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    glBindTexture(GL_TEXTURE_2D, mask1.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (it.x_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    }
    if (it.y_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    glBindTexture(GL_TEXTURE_2D, mask2.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (it.x_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    }
    if (it.y_sym >= 0)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    glBindTexture(GL_TEXTURE_2D, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tex_com.get_shader().use();
    tex_com.get_shader().uniform("x_sym", it.x_sym);
    tex_com.get_shader().uniform("y_sym", it.y_sym);
    tex_com.get_shader().texture("tex1", t1);
    tex_com.get_shader().texture("tex2", t2);
    tex_com.get_shader().texture("mask1", mask1);
    tex_com.get_shader().texture("mask2", mask2);
    tex_com.render();
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.texture, 0);
    glViewport(it.w0 * W, it.h0 * H, w * W, h * H);
    copy.use();
    copy.get_shader().texture("tex", tmp_tex);
    copy.render();
    
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
  }
  glEnable(GL_DEPTH_TEST);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
  return t;
}