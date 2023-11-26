#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/camera.h"
#include "diff_render_interface.h"
#include "tinyEngine/postfx.h"

//we implement the same interface to diff render, but without differentiation
//it can be useful for genetic/memetic algorithms
class NonDiffRender : public IDiffRender
{
public:
  NonDiffRender();
  virtual ~NonDiffRender();
  virtual void init_optimization(const std::vector<Texture> &reference_images, 
                                 Settings render_settings, 
                                 bool save_intermediate_images = false) override;

  virtual float render_and_compare_silhouette(const std::vector<float> &positions, 
                                              const std::vector<CameraSettings> &cameras) override;
  
  virtual const float *get_vertex_grad() const override;
private:
  GLuint fbo, results_buf, vao, vbo;
  Shader render_silhouette, diff_loss_sum;
  PostFx diff_loss;
  int tex_w = 0, tex_h = 0;
  bool save_images = false;
  int iteration = 0;
  std::vector<Texture> reference_textures;
  Texture tex1;
  Texture tex2;
};
