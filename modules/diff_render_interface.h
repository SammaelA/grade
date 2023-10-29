#pragma once
#include <vector>
#include <string>
#include "tinyEngine/camera.h"

class IDiffRender
{
public:
  enum CompareFunction
  {
    MSE
  };
  struct Settings
  {
    int image_w = 128;
    int image_h = 128;
    int samples_per_pixel = 16;
    CompareFunction compare_function = CompareFunction::MSE;
  };

  IDiffRender() = default;
  virtual ~IDiffRender() = default;
  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                 Settings render_settings, 
                                 bool save_intermediate_images = false) = 0;

  virtual float render_and_compare_silhouette(const std::vector<float> &positions, 
                                              const std::vector<CameraSettings> &cameras) = 0;
  
  virtual const float *get_vertex_grad() const = 0;
};
