#pragma once
#include <Python.h>
#include <string>
#include <vector>
#include <array>

class MitsubaInterface
{
public:
  enum MitsubaVariant
  {
    CUDA,
    LLVM
  };

  struct RenderSettings
  {
    RenderSettings(int iw, int ih, int spp, MitsubaVariant mv) : 
                   image_w(iw), image_h(ih), samples_per_pixel(spp), mitsubaVar(mv) {};
    int image_w = 128;
    int image_h = 128;
    int samples_per_pixel = 16;
    MitsubaVariant mitsubaVar = MitsubaVariant::CUDA;
  };
  enum LossFunction
  {
    LOSS_MSE,
    LOSS_MSE_SQRT
  };
  //basic mitsuba initialization, call before any other functions
  void init(const std::string &scripts_dir, const std::string &file_name, RenderSettings render_settings);

  void set_model_max_size(int model_max_size);

  //initialize optimization cycle, set model to compare with. Set loss function and render settings for optimization cycle
  void init_optimization(const std::string &reference_image_dir, LossFunction loss_function, int model_max_size, bool save_intermediate_images);

  //render model and save image to file, for debug purposes
  void render_model_to_file(const std::vector<float> &model, const std::string &image_dir);
  
  //renders model amd compare it with reference set by init_optimization function. Returns loss function value. Saves gradients
  //that are used by compute_final_grad
  float render_and_compare(const std::vector<float> &model);

  //generator_jak size is [FLOAT_PER_VERTEX*params_count*vertex_count], final_grad size is [params_count]
  void compute_final_grad(const std::vector<float> &generator_jac, int params_count, int vertex_count, std::vector<float> &final_grad);

  void finish();
//private:
  void show_errors();
  int get_array_from_ctx_internal(const std::string &name, int buffer_id);//returns loaded array size (in floats)
  void set_array_to_ctx_internal(const std::string &name, int buffer_id, int size);//sends size float from buffer to mitsuba context 
  float render_and_compare_internal();//returns loss function value
  void model_to_ctx(const std::vector<float> &model);
  void clear_buffer(int buffer_id, float val = 0);
  int model_max_size = 0;
  int iteration = 0;
  std::array<float *, 4> buffers = {nullptr, nullptr, nullptr, nullptr};
  PyObject *pModule, *mitsubaContext;
};