#pragma once
#include "diff_render.h"
#include "tinyEngine/postfx.h"

class DiffRenderGPU : public IDiffRender
{
public:
  DiffRenderGPU();
  virtual ~DiffRenderGPU();
  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images = false) override;

  virtual float render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params) override;
  
  virtual const float *get_vertex_grad() const override;
  virtual const float *get_scene_params_grad() const override;

private:
  GLuint fbo, results_buf, edges_buf, pdf_buf, grad_buf;
  Shader render_silhouette, diff_loss_sum, edges_get_pdf, edge_sampling;
  PostFx diff_loss;
  int tex_w = 0, tex_h = 0;
  std::vector<Texture> reference_textures;
  Texture tex1;
  Texture tex2;
  bool save_images = false;
  int iter = 0;

  std::vector<float> pos_grad, cam_grad;
};