#include "gauss_blur_precise.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/graphics_utils.h"
#include <vector>

double gaussian(double sigma, double x)
{
  return (1.0/(sqrt(2*PI)*sigma))*exp(-0.5*SQR(x/sigma));
}

std::vector<float> get_1d_gaussian_filter(float sigma, float precision)
{
  std::vector<float> kernel;
  int max_kernel_size = 16;
  float step = 0.001;
  while (kernel.size() < max_kernel_size)
  {
    double val = 0.0;

    //numerically integrate normal distribution PDF
    for (double x = kernel.size() - 0.5; x<kernel.size() + 0.5; x+=step)
    {
      val += step*gaussian(sigma, x);
    }
    if (val > precision || kernel.empty())
      kernel.push_back(val);
    else
      break;
  }

  //normalize our kernel
  double norm = kernel[0];
  for (int i=1;i<kernel.size();i++)
    norm += 2*kernel[i];
  for (float &val : kernel)
    val/=norm;
  
  return kernel;
}

GaussFilter::GaussFilter(float sigma, float precision):
blur("gaussian_blur_precise.fs")
{
  //check if parameters valid
  assert(sigma>0);
  assert(precision>0);

  //create kernel
  kernel = get_1d_gaussian_filter(sigma, precision);

  //create FBO and SSBO
  fbo = create_framebuffer();
  kernel_buf = create_buffer();

  //bind FBO
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  
  //SSBO always has padding of 16 bytes, so we should use vec4 for every float
  std::vector<glm::vec4> kernel_pad;
  for (float &val : kernel)
    kernel_pad.push_back(glm::vec4(val,0,0,0));
  //fill SSBO and bind it
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, kernel_buf);
  glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec4)*kernel_pad.size(), kernel_pad.data(), GL_STATIC_DRAW);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  //check if FBO valid
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
  
  //unbind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

GaussFilter::~GaussFilter()
{
  delete_buffer(kernel_buf);
  delete_framebuffer(fbo);
}
Texture GaussFilter::perform_gauss_blur(Texture &t)
{
  Texture tex = engine::textureManager->create_texture(t.get_W(), t.get_H());
  gauss_blur(t, tex);
  return tex;
}
void GaussFilter::perform_gauss_blur_inplace(Texture &t)
{
  gauss_blur(t,t);
}
void GaussFilter::gauss_blur(Texture &t, Texture &to)
{
  //check texture type
  assert(t.type == GL_TEXTURE_2D);

  //if we don't have tmp texture or it has wrong size, we need to recreate it
  if (!tmp_tex.is_valid() || tmp_tex.type != t.type || 
      tmp_tex.get_H() != t.get_H() || tmp_tex.get_W() != t.get_W())
  {
    tmp_tex = engine::textureManager->create_texture(t.get_W(), t.get_H());
  }

  //bind FBO and SSBO
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, kernel_buf);

  //first pass from t to tmp_tex
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tmp_tex.texture, 0);
  glViewport(0, 0, t.get_W(), t.get_H());
  //start use postfx
  blur.use();

  //postfx has shader inside. Set textures and uniforms to it. 
  blur.get_shader().texture("tex", t);
  blur.get_shader().uniform("pass",0);
  blur.get_shader().uniform("steps",(int)(kernel.size()-1));
  blur.get_shader().uniform("tex_size",glm::vec2(t.get_W(), t.get_H()));

  //call shader for a full-texture pass
  blur.render();

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  //second pass from tmp_tex to t
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, to.texture, 0);
  glViewport(0, 0, t.get_W(), t.get_H());

  blur.use();
  blur.get_shader().texture("tex", tmp_tex);
  blur.get_shader().uniform("pass",1);
  blur.get_shader().uniform("steps",(int)(kernel.size()-1));
  blur.get_shader().uniform("tex_size",glm::vec2(t.get_W(), t.get_H()));
  blur.render();

  //unbind FBO
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}