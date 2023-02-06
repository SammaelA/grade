#pragma once
#include <string>
#include <vector>
#include <array>
#include "diff_geometry_generation.h"
#include "tinyEngine/camera.h"
struct _object;
typedef _object PyObject;
class MitsubaInterface
{
public:
  enum MitsubaVariant
  {
    CUDA,
    LLVM
  };
  enum RenderStyle
  {
    SILHOUETTE,
    MONOCHROME,
    TEXTURED_CONST
  };
  struct RenderSettings
  {
    RenderSettings() = default;
    RenderSettings(int iw, int ih, int spp, MitsubaVariant mv, RenderStyle rs, std::string _texture_name = "") : 
                   image_w(iw), image_h(ih), samples_per_pixel(spp), mitsubaVar(mv), renderStyle(rs), texture_name(_texture_name) {};
    int image_w = 128;
    int image_h = 128;
    int samples_per_pixel = 16;
    MitsubaVariant mitsubaVar = MitsubaVariant::CUDA;
    RenderStyle renderStyle = RenderStyle::SILHOUETTE;
    std::string texture_name = "white.png";
  };
  enum LossFunction
  {
    LOSS_MSE,
    LOSS_MSE_SQRT,
    LOSS_MIXED
  };

  MitsubaInterface(const std::string &scripts_dir, const std::string &file_name);
  ~MitsubaInterface();

  //basic mitsuba initialization, call before any other functions
  void init_scene_and_settings(RenderSettings render_settings);

  //initialize optimization cycle, set model to compare with. Set loss function and render settings for optimization cycle, includes init_scene_and_settings
  void init_optimization(const std::vector<std::string> &reference_image_dir, LossFunction loss_function, int model_max_size, dgen::ModelLayout opt_ml,
                         RenderSettings render_settings, int cameras_count = 1, bool save_intermediate_images = false);

  //WIP. initialize optimization of texture. includes init_scene_and_settings in it
  void init_optimization_with_tex(const std::vector<std::string> &reference_image_dir, const std::string &initial_texture_name,
                                  LossFunction loss_function, int model_max_size, dgen::ModelLayout opt_ml,
                                  RenderSettings render_settings, int cameras_count = 1, bool save_intermediate_images = false);
  //render model and save image to file, for debug purposes
  void render_model_to_file(const std::vector<float> &model, const std::string &image_dir, const dgen::ModelLayout &ml,
                            const CameraSettings &camera, const std::vector<float> &scene_params);
  
  //renders model amd compare it with reference set by init_optimization function. Returns loss function value. Saves gradients
  //that are used by compute_final_grad
  float render_and_compare(const std::vector<float> &model, const CameraSettings &camera, const std::vector<float> &scene_params,
                           double *timers = nullptr);

  //generator_jak size is [FLOAT_PER_VERTEX*params_count*vertex_count], final_grad size is [params_count]
  void compute_final_grad(const std::vector<float> &generator_jac, int params_count, int vertex_count, std::vector<float> &final_grad);

  void finish();
//private:
  void show_errors();
  void set_model_max_size(int model_max_size);
  void init_optimization_internal(const std::string &function_name, const std::vector<std::string> &reference_image_dir,
                                  LossFunction loss_function, int model_max_size, dgen::ModelLayout opt_ml,
                                  RenderSettings render_settings, int cameras_count, bool save_intermediate_images);
  int get_array_from_ctx_internal(const std::string &name, int buffer_id);//returns loaded array size (in floats)
  void set_array_to_ctx_internal(const std::string &name, int buffer_id, int size);//sends size float from buffer to mitsuba context 
  float render_and_compare_internal(int cameras_count);//returns loss function value
  void model_to_ctx(const std::vector<float> &model, const dgen::ModelLayout &ml);
  void camera_to_ctx(const CameraSettings &camera);
  void clear_buffer(int buffer_id, float val = 0);
  int model_max_size = 0;
  int iteration = 0;
  int cameras_count = 1;
  std::array<float *, 4> buffers = {nullptr, nullptr, nullptr, nullptr};
  std::array<std::string, 4> buffer_names = {"vertex_positions", "vertex_normals", "vertex_texcoords", "camera_params"};
  PyObject *pModule = nullptr, *mitsubaContext = nullptr;
  RenderSettings render_settings;
  dgen::ModelLayout opt_model_layout;
};