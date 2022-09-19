#pragma once
#include <Python.h>
#include <string>
#include <vector>
#include <array>

class MitsubaInterface
{
public:
  void test();
  void init(const std::string &scripts_dir, const std::string &file_name);
  void set_model_max_size(int model_max_size);
  void finish();
//private:
  void show_errors();
  int get_array_from_ctx_internal(const std::string &name, int buffer_id);//returns loaded array size (in floats)
  void set_array_to_ctx_internal(const std::string &name, int buffer_id, int size);//sends size float from buffer to mitsuba context 
  float render_and_compare_internal(int grad_buffer_id);//returns loss function value, saves gradients in buffer with grad_buffer_id
  int model_max_size = 0;
  std::array<float *, 4> buffers = {nullptr, nullptr, nullptr, nullptr};
  PyObject *pModule, *mitsubaContext;
};