#include "mitsuba_python_interaction.h"
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

void MitsubaInterface::init(const std::string &scripts_dir, const std::string &file_name)
{
  //Interpreter initialization
  std::string append_path_str = std::string("sys.path.append(\"")+scripts_dir+"\")";
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("import os");
  PyRun_SimpleString(append_path_str.c_str());
  PyRun_SimpleString("print(sys.path)");
  PyObject *pName;
  pName = PyUnicode_FromString(file_name.c_str());
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule)
    show_errors();
  
  //mitsuba context initialization
  PyObject *initFunc, *initArgs, *basePath;
  basePath = PyUnicode_FromString("resources/mitsuba_data/");
  initArgs = PyTuple_Pack(1, basePath);
  
  initFunc = PyObject_GetAttrString(pModule, (char *)"init");
  if (!initFunc)
    show_errors();
  
  mitsubaContext = PyObject_CallObject(initFunc, initArgs);
  if (!mitsubaContext)
    show_errors();
  
  Py_DECREF(initFunc);
  Py_DECREF(initArgs);
  Py_DECREF(basePath);
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
  //logerr("C++: recieve %d floats from Python", data_floats);
  //for (int i = 0; i < 10; i++)
  //  logerr("%f", (float)(buffers[buffer_id][i]));
  Py_DECREF(args);
  Py_DECREF(func);
  Py_DECREF(params);
  Py_DECREF(params_bytes);
  Py_DECREF(params_name);

  return data_floats;
}

void MitsubaInterface::set_array_to_ctx_internal(const std::string &name, int buffer_id, int size)
{
  //logerr("C++: send %d floats to Python", size);
  //for (int i = 0; i < 10; i++)
  //  logerr("%f", (float)(buffers[buffer_id][i]));

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

float MitsubaInterface::render_and_compare_internal(int grad_buffer_id)
{
  PyObject *pFunc, *pIndex, *pArgs, *pValue;

  pFunc = PyObject_GetAttrString(pModule, (char *)"render");
  if (!pFunc)
    show_errors();
  pIndex = PyLong_FromLong(0);
  pArgs = PyTuple_Pack(2, pIndex, mitsubaContext);
  pValue = PyObject_CallObject(pFunc, pArgs);
  if (!pValue)
    show_errors();
  double result = PyFloat_AsDouble(pValue);
  logerr("val %f", (float)result);

  Py_DECREF(pValue);
  Py_DECREF(pIndex);
  Py_DECREF(pArgs);
  Py_DECREF(pFunc);

  return result;
}

void MitsubaInterface::test()
{
  init("/home/sammael/grade/scripts", "emb_test"); // TODO: replace absolute path
  set_model_max_size(49917);
  int sz = get_array_from_ctx_internal("vertex_positions", 0);
  {
    float *data_initial = buffers[0];
    float *data_to_transport = buffers[1];
    float rotate = 0.0;
    glm::vec3 translate = glm::vec3(0, 0.25, 0);
    glm::mat4 tr_mat = glm::translate(glm::rotate(rotate, glm::vec3(0, 1, 0)), translate);
    for (int j = 0; j < model_max_size; j += 3)
    {
      glm::vec4 res = tr_mat * glm::vec4(data_initial[j], data_initial[j + 1], data_initial[j + 2], 1);
      data_to_transport[j] = res.x;
      data_to_transport[j + 1] = res.y;
      data_to_transport[j + 2] = res.z;
    }
  }
  set_array_to_ctx_internal("vertex_positions", 1, sz);
  float loss = render_and_compare_internal(2);
  finish();
}