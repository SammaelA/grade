#include "diff_render.h"
#include <chrono>
DiffRenderMitsubaDefault::DiffRenderMitsubaDefault(MitsubaInterface &_mi, MitsubaInterface::LossFunction loss_function, 
                                                   MitsubaInterface::ModelInfo model_info, 
                                                   float texture_rec_learing_rate):
  mi(_mi),
  _lf(loss_function),
  _minfo(model_info),
  _tlr(texture_rec_learing_rate)                
{

}
  void DiffRenderMitsubaDefault::init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images)
  {
    if (render_settings.renderStyle == MitsubaInterface::RenderStyle::SILHOUETTE)
      mi.init_optimization(reference_image_dir, _lf, render_settings, _minfo, save_intermediate_images);
    else if (render_settings.renderStyle == MitsubaInterface::RenderStyle::TEXTURED_CONST)
      mi.init_optimization_with_tex(reference_image_dir, _lf, render_settings, _minfo, _tlr, save_intermediate_images);
    else
      logerr("DiffRenderMitsubaDefault supports only SILHOUETTE and TEXTURED_CONST render styles");
  }

  float DiffRenderMitsubaDefault::render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params)
  {
    auto begin = std::chrono::steady_clock::now();
    float res = mi.render_and_compare(model, cameras, scene_params);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Diff. render took = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" <<std::endl;
    return res;
  }

  const float *DiffRenderMitsubaDefault::get_vertex_grad() const
  {
    return mi.get_vertex_grad();
  }
  const float *DiffRenderMitsubaDefault::get_scene_params_grad() const
  {
    return mi.get_camera_params_grad();
  }