#pragma once
#include "upg.h"
#include "tinyEngine/texture.h"
#include "generation.h"
#include "tinyEngine/camera.h"

namespace upg
{
  //structure that contains all data for one given view, like mask and maybe some 
  //depth information in future
  struct ReferenceView
  {
    CameraSettings camera;
    bool fixed_camera = true;
    Texture mask;
    std::vector<Texture> resized_masks;
  };

  //All data that we were able to get about the model that we want to reconstruct
  //before the reconstruction started. If we have synthetic (=create by the generator itself)
  //or 3D model reference, we save it here too to get quality metric after the reconstruction
  struct ReconstructionReference
  {
    ComplexModel model;
    std::vector<ReferenceView> images;
  };

  ReconstructionReference get_reference(const Block &input_blk);
  Texture resize_mask(Texture mask, int tex_w, int tex_h, bool to_binary_mask = false);
  float get_image_based_quality(const ReconstructionReference &reference, const ComplexModel &reconstructed_model);
  float get_model_based_quality(const ReconstructionReference &reference, const ComplexModel &reconstructed_model);
  void  render_model_turntable(const Block &hydra_settings, const ComplexModel &model);
}