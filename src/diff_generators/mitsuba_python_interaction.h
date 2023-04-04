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
    TEXTURED_CONST,
    MONOCHROME_DEMO,
    TEXTURED_DEMO
  };
  struct RenderSettings
  {
    RenderSettings() = default;
    RenderSettings(int iw, int ih, int spp, MitsubaVariant mv, RenderStyle rs, int _cameras_count = 1) : 
                   image_w(iw), image_h(ih), samples_per_pixel(spp), mitsubaVar(mv), renderStyle(rs),
                   cameras_count(_cameras_count) {};
    int image_w = 128;
    int image_h = 128;
    int samples_per_pixel = 16;
    MitsubaVariant mitsubaVar = MitsubaVariant::CUDA;
    RenderStyle renderStyle = RenderStyle::SILHOUETTE;
    int cameras_count = 1;
  };
  struct ModelInfo
  {
    //composite model is made from several parts
    //all parts have the same layout, but different parts
    //have different materials and textures
    struct PartInfo
    {
      std::string name = "main_part";
      std::string texture_name = "white.png";
      std::string material_name = "ceramics";
    };

    dgen::ModelLayout layout;
    std::vector<PartInfo> parts;

    ModelInfo() = default;

    PartInfo *get_part(const std::string &name)
    {
      for (auto &p : parts)
      {
        if (p.name == name)
          return &p;
      }
      return nullptr;
    }

    static ModelInfo simple_mesh(std::string texture_name, std::string material_name)
    {
      ModelInfo mi;
      mi.layout = dgen::ModelLayout();
      mi.parts.push_back(PartInfo{"main_part", texture_name, material_name});
      return mi;
    }
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
  void init_scene_and_settings(RenderSettings render_settings, ModelInfo model_info);

  //initialize optimization cycle, set model to compare with. Set loss function and render settings for optimization cycle, includes init_scene_and_settings
  void init_optimization(const std::vector<std::string> &reference_image_dir, LossFunction loss_function,
                         RenderSettings render_settings, ModelInfo model_info,
                         bool save_intermediate_images = false);

  //WIP. initialize optimization of texture. includes init_scene_and_settings in it
  void init_optimization_with_tex(const std::vector<std::string> &reference_image_dir, LossFunction loss_function, 
                                  RenderSettings render_settings, ModelInfo model_info, 
                                  float texture_rec_learing_rate = 0.25,
                                  bool save_intermediate_images = false);
  //render model and save image to file, for debug purposes
  void render_model_to_file(const dgen::DFModel &model, const std::string &image_dir,
                            const CameraSettings &camera, const std::vector<float> &scene_params);
  
  //renders model amd compare it with reference set by init_optimization function. Returns loss function value. Saves gradients
  //that are used by compute_final_grad
  float render_and_compare(const dgen::DFModel &model, const CameraSettings &camera, const std::vector<float> &scene_params,
                           double *timers = nullptr);

  //render model from different angles, merge them into one image and save it to file, for debug purposes
  //scene_params SHOULD NOT rotate or translate the model
  void render_multicam_demo(RenderSettings render_settings, ModelInfo model_info, 
                            const dgen::DFModel &model, const std::string &image_dir,
                            const std::vector<float> &scene_params, const CameraSettings &camera,
                            int rotations_x = 4, int rotations_y = 1);

  //generator_jak size is [FLOAT_PER_VERTEX*params_count*vertex_count], final_grad size is [params_count]
  void compute_final_grad(const std::vector<float> &generator_jac, int params_count, int vertex_count, std::vector<float> &final_grad);

  void finish();

  static std::vector<std::string> get_all_available_materials();
  static std::string get_default_material();
  static CameraSettings get_camera_from_scene_params(const std::vector<float> &scene_params);
//private:
  void show_errors();
  void set_model_max_size(int model_max_size);
  void init_optimization_internal(const std::string &function_name, const std::vector<std::string> &reference_image_dir,
                                  LossFunction loss_function, RenderSettings render_settings, ModelInfo model_info,
                                  float texture_rec_learing_rate, bool save_intermediate_images);
  int get_array_from_ctx_internal(const std::string &name, int buffer_id);//returns loaded array size (in floats)
  void set_array_to_ctx_internal(const std::string &name, int buffer_id, int size);//sends size float from buffer to mitsuba context 
  float render_and_compare_internal();//returns loss function value
  void model_to_ctx(const dgen::DFModel &model);
  void camera_to_ctx(const CameraSettings &camera);
  void clear_buffer(int buffer_id, float val = 0);
  int get_camera_buffer_id()
  {
    return buffers.size() - 1;
  }
  int model_max_size = 0;
  int iteration = 0;
  std::vector<float *> buffers;
  std::vector<std::string> buffer_names;
  std::vector<int> active_parts;
  PyObject *pModule = nullptr, *mitsubaContext = nullptr;
  RenderSettings render_settings;
  ModelInfo model_info;
};