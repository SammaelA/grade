#include "mitsuba_python_interaction.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/silhouette.h"
#include "diff_geometry_generation.h"
#include <iostream>
#include "common_utils/utility.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>

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

  Py_DECREF(mitsubaContext);
  Py_DECREF(pModule);
  Py_Finalize();
}

void MitsubaInterface::init(const std::string &scripts_dir, const std::string &file_name, RenderSettings render_settings)
{
  //Interpreter initialization
  std::string append_path_str = std::string("sys.path.append(\"")+scripts_dir+"\")";
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("import os");
  PyRun_SimpleString(append_path_str.c_str());
  PyObject *pName;
  pName = PyUnicode_FromString(file_name.c_str());
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule)
    show_errors();
  
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
  PyObject *initFunc, *initArgs, *basePath, *iw_arg, *ih_arg, *spp_arg, *mv;
  basePath = PyUnicode_FromString("resources/mitsuba_data/");
  iw_arg = PyLong_FromLong(render_settings.image_w);
  ih_arg = PyLong_FromLong(render_settings.image_h);
  spp_arg = PyLong_FromLong(render_settings.samples_per_pixel);
  mv = PyUnicode_FromString(mitsuba_var.c_str());
  initArgs = PyTuple_Pack(5, basePath, iw_arg, ih_arg, spp_arg, mv);
  
  initFunc = PyObject_GetAttrString(pModule, (char *)"init");
  if (!initFunc)
    show_errors();
  
  mitsubaContext = PyObject_CallObject(initFunc, initArgs);
  if (!mitsubaContext)
    show_errors();
  
  Py_DECREF(initFunc);
  Py_DECREF(initArgs);
  Py_DECREF(basePath);
  Py_DECREF(iw_arg);
  Py_DECREF(ih_arg);
  Py_DECREF(spp_arg);
  Py_DECREF(mv);
}

void MitsubaInterface::init_optimization(const std::string &reference_image_dir, LossFunction loss_function, int model_max_size, bool save_intermediate_images)
{
  std::string loss_function_name = "F_loss";
  switch (loss_function)
  {
  case LossFunction::LOSS_MSE :
    loss_function_name = "F_loss";
    break;

  case LossFunction::LOSS_MSE_SQRT :
    loss_function_name = "F_loss_sqrt";
    break;

  default:
    loss_function_name = "F_loss";
    break;
  }
  PyObject *func, *args, *ref_dir_arg, *func_ret, *loss_func, *int_im;

  func = PyObject_GetAttrString(pModule, (char *)"init_optimization");
  ref_dir_arg = PyUnicode_FromString(reference_image_dir.c_str());
  loss_func = PyObject_GetAttrString(pModule, loss_function_name.c_str());
  int_im = PyLong_FromLong((int)save_intermediate_images);
  args = PyTuple_Pack(4, mitsubaContext, ref_dir_arg, loss_func, int_im);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  set_model_max_size(model_max_size);
  iteration = 0;

  Py_DECREF(func);
  Py_DECREF(args);
  Py_DECREF(ref_dir_arg);
  Py_DECREF(func_ret);
  Py_DECREF(loss_func);
  Py_DECREF(int_im);
}

void MitsubaInterface::model_to_ctx(const std::vector<float> &model)
{
  int vertex_count = model.size() / FLOAT_PER_VERTEX;
  clear_buffer(0, 0.0f);
  clear_buffer(1, 1.0f);
  clear_buffer(2, 0.0f);
  for (int i = 0; i < vertex_count; i++)
  {
    buffers[0][3 * i] = model[FLOAT_PER_VERTEX * i];
    buffers[0][3 * i + 1] = model[FLOAT_PER_VERTEX * i + 1];
    buffers[0][3 * i + 2] = model[FLOAT_PER_VERTEX * i + 2];

    buffers[1][3 * i] = model[FLOAT_PER_VERTEX * i + 3];
    buffers[1][3 * i + 1] = model[FLOAT_PER_VERTEX * i + 4];
    buffers[1][3 * i + 2] = model[FLOAT_PER_VERTEX * i + 5];

    buffers[2][2 * i] = model[FLOAT_PER_VERTEX * i + 6];
    buffers[2][2 * i + 1] = model[FLOAT_PER_VERTEX * i + 7];
  }

  set_array_to_ctx_internal("vertex_positions", 0, 3 * vertex_count);
  set_array_to_ctx_internal("vertex_normals", 1, 3 * vertex_count);
  set_array_to_ctx_internal("vertex_texcoords", 2, 2 * vertex_count);
}

void MitsubaInterface::render_model_to_file(const std::vector<float> &model, const std::string &image_dir)
{
  model_to_ctx(model);

  PyObject *func, *args, *ref_dir_arg, *func_ret;

  func = PyObject_GetAttrString(pModule, (char *)"render_and_save_to_file");
  ref_dir_arg = PyUnicode_FromString(image_dir.c_str());
  args = PyTuple_Pack(2, mitsubaContext, ref_dir_arg);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  Py_DECREF(func);
  Py_DECREF(args);
  Py_DECREF(ref_dir_arg);
  Py_DECREF(func_ret);

  Texture t = engine::textureManager->load_unnamed_tex(image_dir);
  SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
  Texture tex = se.get_silhouette(t);
  engine::textureManager->save_png_directly(tex, image_dir);
}

float MitsubaInterface::render_and_compare(const std::vector<float> &model)
{
  model_to_ctx(model);
  float loss = render_and_compare_internal();
  get_array_from_ctx_internal("vertex_positions_grad", 0);
  get_array_from_ctx_internal("vertex_normals_grad", 1);
  get_array_from_ctx_internal("vertex_texcoords_grad", 2);
  return loss;
}

void MitsubaInterface::compute_final_grad(const std::vector<float> &jac, int params_count, int vertex_count, 
                                          std::vector<float> &final_grad)
{
  for (int i = 0; i < vertex_count; i++)
  {
    for (int j = 0; j < params_count; j++)
    {
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i) * params_count + j] * buffers[0][3 * i];
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 1) * params_count + j] * buffers[0][3 * i + 1];
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 2) * params_count + j] * buffers[0][3 * i + 2];

      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 3) * params_count + j] * buffers[1][3 * i];
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 4) * params_count + j] * buffers[1][3 * i + 1];
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 5) * params_count + j] * buffers[1][3 * i + 2];

      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 6) * params_count + j] * buffers[2][2 * i];
      final_grad[j] += jac[(FLOAT_PER_VERTEX * i + 7) * params_count + j] * buffers[2][2 * i + 1];
    }
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
      buffers[i] = new float[model_max_size];
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
  Py_DECREF(args);
  Py_DECREF(func);
  Py_DECREF(params);
  Py_DECREF(params_bytes);
  Py_DECREF(params_name);

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

  Py_DECREF(args);
  Py_DECREF(func);
  Py_DECREF(params);
  Py_DECREF(params_n);
  Py_DECREF(params_bytes);
  Py_DECREF(params_name);
}

float MitsubaInterface::render_and_compare_internal()
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

  Py_DECREF(pValue);
  Py_DECREF(pIndex);
  Py_DECREF(pArgs);
  Py_DECREF(pFunc);

  return result;
}

void MitsubaInterface::clear_buffer(int buffer_id, float val)
{
  std::fill_n(buffers[buffer_id], model_max_size, val);
}