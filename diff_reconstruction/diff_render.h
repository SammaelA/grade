#pragma once
#include "mitsuba_python_interaction.h"

class IDiffRender
{
public:
  IDiffRender() = default;
  virtual ~IDiffRender() = default;
  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images = false) = 0;

  virtual float render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params) = 0;
  
  virtual const float *get_vertex_grad() const = 0;
  virtual const float *get_scene_params_grad() const = 0;
};

class DiffRenderMitsubaDefault : public IDiffRender
{
public:
  DiffRenderMitsubaDefault(MitsubaInterface &mi, MitsubaInterface::LossFunction loss_function, 
                           MitsubaInterface::ModelInfo model_info, 
                           float texture_rec_learing_rate);
  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images = false) override;

  virtual float render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params) override;
  
  virtual const float *get_vertex_grad() const override;
  virtual const float *get_scene_params_grad() const override;

private:
  MitsubaInterface &mi;
  MitsubaInterface::LossFunction _lf;
  MitsubaInterface::ModelInfo _minfo;
  float _tlr;

};