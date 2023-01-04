#include <Python.h>
#include "common_utils/python_engine.h"
#include "mitsuba_python_interaction.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/silhouette.h"
#include "diff_geometry_generation.h"
#include <iostream>
#include "common_utils/utility.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>

#define DEL(X) if (X) {Py_DECREF(X);}
void MitsubaInterface::show_errors()
{
  PyObject *pExcType, *pExcValue, *pExcTraceback;
  PyErr_Fetch(&pExcType, &pExcValue, &pExcTraceback);
  if (pExcType != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcType);
    logerr("An error occurred:");
    logerr("- EXC type: %s", PyUnicode_AsUTF8(pRepr));
    Py_DecRef(pRepr);
    Py_DecRef(pExcType);
  }
  if (pExcValue != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcValue);
    logerr("An error occurred:");
    logerr("- EXC value: %s", PyUnicode_AsUTF8(pRepr));
    Py_DecRef(pRepr);
    Py_DecRef(pExcValue);
  }
  if (pExcTraceback != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcValue);
    logerr("An error occurred:");
    logerr("- EXC traceback: %s", PyUnicode_AsUTF8(pRepr));
    Py_DecRef(pRepr);
    Py_DecRef(pExcTraceback);
  }
}

void MitsubaInterface::finish()
{
  for (int i = 0; i < buffers.size(); i++)
  {
    if (buffers[i])
    {
      delete[] buffers[i];
      buffers[i] = nullptr;
    }
  }
}

MitsubaInterface::~MitsubaInterface()
{
  DEL(mitsubaContext);
  DEL(pModule);
  Py_Finalize();
}

MitsubaInterface::MitsubaInterface(const std::string &scripts_dir, const std::string &file_name)
{
  //Interpreter initialization
  std::string append_path_str = std::string("sys.path.append(\"")+scripts_dir+"\")";
  python_engine::init();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("import os");
  PyRun_SimpleString(append_path_str.c_str());
  PyObject *pName;
  pName = PyUnicode_FromString(file_name.c_str());
  pModule = PyImport_Import(pName);
  DEL(pName);
  if (!pModule)
    show_errors();
}

void MitsubaInterface::init_scene_and_settings(RenderSettings _render_settings)
{
  render_settings = _render_settings;
  //mitsuba context initialization
  std::string mitsuba_var = "";
  switch (render_settings.mitsubaVar)
  {
  case MitsubaVariant::CUDA :
    mitsuba_var = "cuda_ad_rgb";
    break;
  case MitsubaVariant::LLVM :
    mitsuba_var = "llvm_ad_rgb";
    break;
  default:
    mitsuba_var = "cuda_ad_rgb";
    break;
  }

  std::string render_style = "";
  switch (render_settings.renderStyle)
  {
  case RenderStyle::SILHOUETTE:
    render_style = "silhouette";
    break;
  case RenderStyle::MONOCHROME:
    render_style = "monochrome";
    break;
  case RenderStyle::TEXTURED_CONST:
    render_style = "textured_const";
    break;
  default:
    render_style = "silhouette";
    break;
  }
  PyObject *initFunc, *initArgs, *basePath, *iw_arg, *ih_arg, *spp_arg, *mv, *rs, *tn;
  basePath = PyUnicode_FromString("resources/mitsuba_data/");
  iw_arg = PyLong_FromLong(render_settings.image_w);
  ih_arg = PyLong_FromLong(render_settings.image_h);
  spp_arg = PyLong_FromLong(render_settings.samples_per_pixel);
  mv = PyUnicode_FromString(mitsuba_var.c_str());
  rs = PyUnicode_FromString(render_style.c_str());
  tn = PyUnicode_FromString(render_settings.texture_name.c_str());
  initArgs = PyTuple_Pack(7, basePath, iw_arg, ih_arg, spp_arg, mv, rs, tn);
  
  initFunc = PyObject_GetAttrString(pModule, (char *)"init");
  if (!initFunc)
    show_errors();
  
  if (mitsubaContext)
    DEL(mitsubaContext);

  mitsubaContext = PyObject_CallObject(initFunc, initArgs);
  if (!mitsubaContext)
    show_errors();
  
  DEL(initFunc);
  DEL(initArgs);
  DEL(basePath);
  DEL(iw_arg);
  DEL(ih_arg);
  DEL(spp_arg);
  DEL(mv);
  DEL(rs);
}

std::string get_loss_function_name(MitsubaInterface::LossFunction loss_function)
{
  std::string loss_function_name = "F_loss_mse";
  switch (loss_function)
  {
  case MitsubaInterface::LossFunction::LOSS_MSE :
    loss_function_name = "F_loss_mse";
    break;

  case MitsubaInterface::LossFunction::LOSS_MSE_SQRT :
    loss_function_name = "F_loss_mse_sqrt";
    break;
  
  case MitsubaInterface::LossFunction::LOSS_MIXED :
    loss_function_name = "F_loss_mixed";
    break;

  default:
    loss_function_name = "F_loss_mse";
    break;
  }
  return loss_function_name;
}
void MitsubaInterface::init_optimization(const std::vector<std::string> &reference_image_dir, LossFunction loss_function, int model_max_size, 
                                         dgen::ModelLayout opt_ml,
                                         RenderSettings render_settings, int cam_count, bool save_intermediate_images)
{
  init_optimization_internal("init_optimization", reference_image_dir, loss_function, model_max_size, opt_ml, 
                             render_settings, cam_count, save_intermediate_images);
}

void MitsubaInterface::init_optimization_internal(const std::string &function_name,
                                                  const std::vector<std::string> &reference_images_dir, LossFunction loss_function, int model_max_size, 
                                                  dgen::ModelLayout opt_ml, RenderSettings render_settings, 
                                                  int cam_count, bool save_intermediate_images)
{
  init_scene_and_settings(render_settings);

  opt_model_layout = opt_ml;
  cameras_count = cam_count;
  std::string loss_function_name = get_loss_function_name(loss_function);
  
  //save all strings as "ref1.png#ref2.png#ref3.png"
  std::string full_ref_string = "";
  assert(reference_images_dir.size() >= cam_count);
  for (int i=0;i<cam_count;i++)
  {
    full_ref_string += reference_images_dir[i];
    if (i < cam_count - 1)
     full_ref_string += "#"; 
  }

  PyObject *func, *args, *ref_dir_arg, *func_ret, *loss_func, *c_cnt, *int_im;

  func = PyObject_GetAttrString(pModule, function_name.c_str());
  ref_dir_arg = PyUnicode_FromString(full_ref_string.c_str());
  loss_func = PyObject_GetAttrString(pModule, loss_function_name.c_str());
  int_im = PyLong_FromLong((int)save_intermediate_images);
  c_cnt = PyLong_FromLong(cam_count);
  args = PyTuple_Pack(5, mitsubaContext, ref_dir_arg, loss_func, c_cnt, int_im);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  set_model_max_size(model_max_size);
  iteration = 0;

  DEL(func);
  DEL(args);
  DEL(ref_dir_arg);
  DEL(func_ret);
  DEL(loss_func);
  DEL(c_cnt);
  DEL(int_im);
}

void MitsubaInterface::init_optimization_with_tex(const std::vector<std::string> &reference_image_dir, const std::string &initial_texture_name,
                                                  LossFunction loss_function, int model_max_size, dgen::ModelLayout opt_ml,
                                                  RenderSettings render_settings, int cam_count, bool save_intermediate_images)
{
  render_settings.renderStyle = RenderStyle::TEXTURED_CONST;
  render_settings.texture_name = initial_texture_name;
  
  init_optimization_internal("init_optimization_with_tex", reference_image_dir, loss_function, model_max_size, opt_ml, 
                             render_settings, cam_count, save_intermediate_images);
}

void MitsubaInterface::model_to_ctx(const std::vector<float> &model, const dgen::ModelLayout &ml)
{
  int vertex_count = model.size() / ml.f_per_vert;
  assert(ml.offsets.size() - 1 <= buffers.size());
  for (int i=0;i<ml.offsets.size() - 1;i++)
  {
    int offset = ml.offsets[i];
    int size = ml.offsets[i + 1] - ml.offsets[i];
    if (offset >= 0 && size > 0)
    {
      clear_buffer(i, 0.0f);
      for (int j = 0; j < vertex_count; j++)
        memcpy(buffers[i] + size*j, model.data() + ml.f_per_vert * j + offset, sizeof(float)*size);
      set_array_to_ctx_internal(buffer_names[i], i, size * vertex_count);
    }
  }
  show_errors();
}

void MitsubaInterface::camera_to_ctx(const CameraSettings &camera)
{
  PyObject *func, *args, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *p9, *p10, *p11, *p12, *func_ret;

  func = PyObject_GetAttrString(pModule, (char *)"set_camera");

  p1 = PyFloat_FromDouble(camera.origin.x);
  p2 = PyFloat_FromDouble(camera.origin.y);
  p3 = PyFloat_FromDouble(camera.origin.z);

  p4 = PyFloat_FromDouble(camera.target.x);
  p5 = PyFloat_FromDouble(camera.target.y);
  p6 = PyFloat_FromDouble(camera.target.z);

  p7 = PyFloat_FromDouble(camera.up.x);
  p8 = PyFloat_FromDouble(camera.up.y);
  p9 = PyFloat_FromDouble(camera.up.z);

  p10 = PyFloat_FromDouble(180 * camera.fov_rad / PI);
  p11 = PyLong_FromLong(render_settings.image_w);
  p12 = PyLong_FromLong(render_settings.image_h);

  args = PyTuple_Pack(13, mitsubaContext, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  DEL(func);
  DEL(args);
  DEL(p1);
  DEL(p2);
  DEL(p3);
  DEL(p4);
  DEL(p5);
  DEL(p6);
  DEL(p7);
  DEL(p8);
  DEL(p9);
  DEL(p10);
  DEL(p11);
  DEL(p12);
  DEL(func_ret);
}

void MitsubaInterface::render_model_to_file(const std::vector<float> &model, const std::string &image_dir, const dgen::ModelLayout &ml,
                                            const CameraSettings &camera)
{
  if (model_max_size < model.size()/ml.f_per_vert)
    set_model_max_size(model.size()/ml.f_per_vert);
  
  model_to_ctx(model, ml);
  camera_to_ctx(camera);

  PyObject *func, *args, *ref_dir_arg, *func_ret;

  func = PyObject_GetAttrString(pModule, (char *)"render_and_save_to_file");
  ref_dir_arg = PyUnicode_FromString(image_dir.c_str());
  args = PyTuple_Pack(2, mitsubaContext, ref_dir_arg);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  DEL(func);
  DEL(args);
  DEL(ref_dir_arg);
  DEL(func_ret);
}

float MitsubaInterface::render_and_compare(const std::vector<float> &model, const CameraSettings &camera, const std::vector<float> &camera_params,
                                           double *timers)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  model_to_ctx(model, opt_model_layout);
  camera_to_ctx(camera);

  //camera params contains 6*n floats: (pos.x, pos.y, pos.z, rot_x, rot_y, rot_z) for each camera
  assert(camera_params.size() == 6 * cameras_count);
  constexpr int cameras_buf_n = 3;

  if (camera_params.empty())
    std::fill_n(buffers[cameras_buf_n], 6*cameras_count, 0);
  else
    memcpy(buffers[3], camera_params.data(), sizeof(float)*6*cameras_count);
  set_array_to_ctx_internal(buffer_names[cameras_buf_n], cameras_buf_n, 6*cameras_count);
  show_errors();

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  float loss = render_and_compare_internal(cameras_count);
  std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
  for (int i=0;i<opt_model_layout.offsets.size() - 1;i++)
  {
    int offset = opt_model_layout.offsets[i];
    int size = opt_model_layout.offsets[i + 1] - opt_model_layout.offsets[i];
    if (offset >= 0 && size > 0)
      get_array_from_ctx_internal(buffer_names[i] + "_grad", i);
  }
  get_array_from_ctx_internal(buffer_names[cameras_buf_n] + "_grad", cameras_buf_n);
  std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
  if (timers)
  {
    timers[2] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    timers[3] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    timers[4] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
  }
  return loss;
}

void MitsubaInterface::compute_final_grad(const std::vector<float> &jac, int params_count, int vertex_count, 
                                          std::vector<float> &final_grad)
{
  for (int off=0;off<opt_model_layout.offsets.size() - 1;off++)
  {
    int offset = opt_model_layout.offsets[off];
    int size = opt_model_layout.offsets[off + 1] - opt_model_layout.offsets[off];
    if (offset >= 0 && size > 0)
    {
      for (int i = 0; i < vertex_count; i++)
      {
        for (int j = 0; j < params_count; j++)
        {
          for (int k = 0; k < size; k++)
          {
            final_grad[j] += jac[(opt_model_layout.f_per_vert * i + offset + k) * params_count + j] * buffers[off][size * i + k];
          }
        }
      }
    }

    for (int i=params_count; i< final_grad.size(); i++)
      final_grad[i] = buffers[3][i - params_count]; // gradient by camera params
  }
}

void MitsubaInterface::set_model_max_size(int _model_max_size)
{
  model_max_size = _model_max_size;
  if (model_max_size >= 0)
  {
    for (int i=0;i<buffers.size(); i++)
    {
      if (buffers[i])
        delete[] buffers[i];
      buffers[i] = new float[4*model_max_size];
    }
  }
}

int MitsubaInterface::get_array_from_ctx_internal(const std::string &name, int buffer_id)
{
  PyObject *func, *args, *params, *params_bytes, *params_name;
  params_name = PyUnicode_FromString(name.c_str());
  args = PyTuple_Pack(2, mitsubaContext, params_name);
  func = PyObject_GetAttrString(pModule, (char *)"get_params");
  params = PyObject_CallObject(func, args);
  if (!params)
    show_errors();
  params_bytes = PyObject_Bytes(params);
  if (!params_bytes)
    show_errors();

  int sz = PyBytes_Size(params_bytes);
  int data_floats = sz / sizeof(float);
  if (data_floats > model_max_size)
  {
    logerr("Python array %s contains %d float, while buffer size is %d. Some data will be ignored", name.c_str(), data_floats, model_max_size);
  }
  char *data = PyBytes_AsString(params_bytes);
  memcpy(buffers[buffer_id], data, MIN(sz, model_max_size*sizeof(float)));
  DEL(args);
  DEL(func);
  DEL(params);
  DEL(params_bytes);
  DEL(params_name);

  return data_floats;
}

void MitsubaInterface::set_array_to_ctx_internal(const std::string &name, int buffer_id, int size)
{
  PyObject *func, *args, *params_n, *params_bytes, *params, *params_name;
  params_name = PyUnicode_FromString(name.c_str());
  params_n = PyLong_FromLong(size);
  params_bytes = PyBytes_FromStringAndSize((const char *)buffers[buffer_id], sizeof(float) * size);
  args = PyTuple_Pack(4, mitsubaContext, params_name, params_bytes, params_n);
  func = PyObject_GetAttrString(pModule, (char *)"set_params");
  params = PyObject_CallObject(func, args);

  DEL(args);
  DEL(func);
  DEL(params);
  DEL(params_n);
  DEL(params_bytes);
  DEL(params_name);
}

float MitsubaInterface::render_and_compare_internal(int cameras_count)
{
  PyObject *pFunc, *pIndex, *pArgs, *pValue;

  pFunc = PyObject_GetAttrString(pModule, (char *)"render");
  if (!pFunc)
    show_errors();
  pIndex = PyLong_FromLong(iteration);
  iteration++;
  pArgs = PyTuple_Pack(2, pIndex, mitsubaContext);
  pValue = PyObject_CallObject(pFunc, pArgs);
  if (!pValue)
    show_errors();
  double result = PyFloat_AsDouble(pValue);

  DEL(pValue);
  DEL(pIndex);
  DEL(pArgs);
  DEL(pFunc);

  return result;
}

void MitsubaInterface::clear_buffer(int buffer_id, float val)
{
  std::fill_n(buffers[buffer_id], model_max_size, val);
}