#include "depth_extract_compare.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include <glm/glm.hpp>
#include <vector>

DepthLossCalculator::DepthLossCalculator():
render_model({"uv_coords.vs", "uv_coords.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
depth_postprocess("depth_postprocess.fs")
{
  fbo = create_framebuffer();
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

DepthLossCalculator::~DepthLossCalculator()
{
  delete_framebuffer(fbo);
}

Texture DepthLossCalculator::get_depth(Model &m, const CameraSettings &camera, int w, int h)
{
  Texture color = engine::textureManager->create_texture(w, h, GL_RGBA32F, 1, NULL, GL_RGBA, GL_FLOAT);
  Texture depth = engine::textureManager->create_texture(w, h, GL_DEPTH_COMPONENT32, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);

  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth.texture, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color.texture, 0);
  glViewport(0, 0, w, h);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glm::mat4 projection = glm::perspective(camera.fov_rad, (float)w/h, camera.z_near, camera.z_far);
  glm::mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);
  render_model.use();

  render_model.uniform("projection", projection);
  render_model.uniform("view", view);
  m.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color.texture, 0);
  glViewport(0, 0, w, h);
  depth_postprocess.use();
  depth_postprocess.get_shader().texture("depth_tex", depth);
  depth_postprocess.get_shader().uniform("tex_size", glm::vec2(w, h));
  depth_postprocess.get_shader().uniform("zNear", camera.z_near);
  depth_postprocess.get_shader().uniform("zFar", camera.z_far);
  depth_postprocess.render();
  
  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

  return color;
}

float DepthLossCalculator::get_loss(Model &m, Texture reference_depth, const CameraSettings &camera)
{
  Texture d = get_depth(m, camera, reference_depth.get_W(), reference_depth.get_H());

  float *data_1 = new float[4*d.get_W()*d.get_H()];
  float *data_2 = new float[4*d.get_W()*d.get_H()];

  glBindTexture(GL_TEXTURE_2D, d.texture);
  glGetTexImage(GL_TEXTURE_2D,
                0,
                GL_RGBA,
                GL_FLOAT,
                data_1);
  glBindTexture(GL_TEXTURE_2D, reference_depth.texture);
  glGetTexImage(GL_TEXTURE_2D,
                0,
                GL_RGBA,
                GL_FLOAT,
                data_2);
  glBindTexture(d.type, 0);

  double diff = 0.0;
  double sum = 1e-6;
  for (int i=0;i<4*d.get_W()*d.get_H();i+=4)
  {
    if (data_1[i+3] > 0.5 || data_2[i+3] > 0.5)
    {
      diff += abs(data_1[i] - data_2[i]); //R channel where silhouettes overlap
      sum += 1;
      //debug("%2d ", MIN((int)(1000*abs(data_1[i] - data_2[i])), 99));
    }
    //else
    //  debug("00 ");
    //if ((i+4) % (4*d.get_W()) == 0)
    //  debugnl();
  }
  //debugnl();
  delete[] data_1;
  delete[] data_2;

  return diff / sum;
}