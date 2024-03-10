#include <Python.h>
#include "python_engine.h"
#include "mitsuba_python_interaction.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/silhouette.h"
#include "diff_geometry_generation.h"
#include <iostream>
#include "common_utils/utility.h"
#include "common_utils/matrix_transform.h"
#include "tinyEngine/engine.h"

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
    PyObject *pRepr = PyObject_Repr(pExcTraceback);
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
  buffers.clear();
  buffer_names.clear();
  model_max_size = 0;
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

void MitsubaInterface::init_scene_and_settings(RenderSettings _render_settings, ModelInfo _model_info)
{
  finish();
  render_settings = _render_settings;
  model_info = _model_info;

  int mesh_parts_count = model_info.parts.size();
  for (int part_id = 0; part_id < mesh_parts_count; part_id++)
  {
    buffer_names.push_back("vertex_positions_"+std::to_string(part_id));
    buffers.push_back(nullptr);

    buffer_names.push_back("vertex_normals_"+std::to_string(part_id));
    buffers.push_back(nullptr);

    buffer_names.push_back("vertex_texcoords_"+std::to_string(part_id));
    buffers.push_back(nullptr);
  }
  buffer_names.push_back("camera_params");
  buffers.push_back(nullptr);

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
  case RenderStyle::TEXTURED_DEMO:
    render_style = "textured_demo";
    break;
  case RenderStyle::MONOCHROME_DEMO:
    render_style = "monochrome_demo";
    break;
  default:
    render_style = "silhouette";
    break;
  }

  std::string texture_names = "";
  std::string material_names = "";
  for (int i=0;i<model_info.parts.size();i++)
  {
    if (i != 0)
    {
      texture_names +="|";
      material_names +="|";
    }
    texture_names += model_info.parts[i].texture_name;
    material_names += model_info.parts[i].material_name;
  }

  PyObject *initFunc, *initArgs, *basePath, *iw_arg, *ih_arg, *spp_arg, *mv, *rs, *tn, *mn;
  basePath = PyUnicode_FromString("resources/mitsuba_data/");
  iw_arg = PyLong_FromLong(render_settings.image_w);
  ih_arg = PyLong_FromLong(render_settings.image_h);
  spp_arg = PyLong_FromLong(render_settings.samples_per_pixel);
  mv = PyUnicode_FromString(mitsuba_var.c_str());
  rs = PyUnicode_FromString(render_style.c_str());
  tn = PyUnicode_FromString(texture_names.c_str());
  mn = PyUnicode_FromString(material_names.c_str());
  initArgs = PyTuple_Pack(8, basePath, iw_arg, ih_arg, spp_arg, mv, rs, tn, mn);
  
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
  DEL(tn);
  DEL(mn);
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
void MitsubaInterface::init_optimization(const std::vector<std::string> &reference_image_dir, LossFunction loss_function,
                                         RenderSettings render_settings, ModelInfo model_info,
                                         bool save_intermediate_images)
{
  init_optimization_internal("init_optimization", reference_image_dir, loss_function, render_settings, model_info,
                             0, save_intermediate_images);
}

void MitsubaInterface::init_optimization_internal(const std::string &function_name, const std::vector<std::string> &reference_images_dir,
                                                  LossFunction loss_function, RenderSettings render_settings, ModelInfo model_info,
                                                  float texture_rec_learing_rate, bool save_intermediate_images)
{
  init_scene_and_settings(render_settings, model_info);

  std::string loss_function_name = get_loss_function_name(loss_function);
  
  //save all strings as "ref1.png#ref2.png#ref3.png"
  std::string full_ref_string = "";
  for (int i=0;i<reference_images_dir.size();i++)
  {
    full_ref_string += reference_images_dir[i];
    if (i < reference_images_dir.size() - 1)
     full_ref_string += "#"; 
  }

  PyObject *func, *args, *ref_dir_arg, *func_ret, *loss_func, *lr, *c_cnt, *int_im;

  func = PyObject_GetAttrString(pModule, function_name.c_str());
  ref_dir_arg = PyUnicode_FromString(full_ref_string.c_str());
  loss_func = PyObject_GetAttrString(pModule, loss_function_name.c_str());
  lr = PyFloat_FromDouble(texture_rec_learing_rate);
  int_im = PyLong_FromLong((int)save_intermediate_images);
  c_cnt = PyLong_FromLong(reference_images_dir.size());
  args = PyTuple_Pack(6, mitsubaContext, ref_dir_arg, loss_func, lr, c_cnt, int_im);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  iteration = 0;

  DEL(func);
  DEL(args);
  DEL(ref_dir_arg);
  DEL(func_ret);
  DEL(loss_func);
  DEL(lr);
  DEL(c_cnt);
  DEL(int_im);
}

void MitsubaInterface::init_optimization_with_tex(const std::vector<std::string> &reference_image_dir, LossFunction loss_function, 
                                                  RenderSettings render_settings, ModelInfo model_info, 
                                                  float texture_rec_learing_rate,
                                                  bool save_intermediate_images)
{
  render_settings.renderStyle = RenderStyle::TEXTURED_CONST;
  init_optimization_internal("init_optimization_with_tex", reference_image_dir, loss_function, render_settings,
                             model_info, texture_rec_learing_rate, save_intermediate_images);
}

void MitsubaInterface::model_to_ctx(const dgen::DFModel &model)
{
  active_parts.clear();
  int start_buffer_offset = 0;
  dgen::PartOffsets off = model.second;
  off.push_back({"", model.first.size()});//to make offset calculations simplier

  for (auto &part : model_info.parts)
  {
    int part_size = 0;
    const float *part_data = nullptr;
    for (int i=0;i<off.size()-1;i++)
    {
      if (off[i].first == part.name)
      {
        part_data = model.first.data() + off[i].second;
        part_size = off[i+1].second - off[i].second;
      }
    }

    float placeholder_triangle[] = 
    {
      0,0,0, 0,0,1, 0,0,
      0.001,0,0, 0,0,1, 0,1,
      0,0.001,0, 0,0,1, 1,0
    };
    if (part_data == nullptr || part_size <= 0)
    {
      //This part does not exist in model (it's OK) or corrupted
      part_data = placeholder_triangle;
      part_size = sizeof(placeholder_triangle)/sizeof(float);
    }
    else
    {
      active_parts.push_back(start_buffer_offset/3);
    }

    auto &ml = model_info.layout;
    int vertex_count = part_size / ml.f_per_vert;
    if (model_max_size < vertex_count)
      set_model_max_size(vertex_count);
    assert(start_buffer_offset + ml.offsets.size() - 1 <= buffers.size());
    for (int i=0;i<ml.offsets.size() - 1;i++)
    {
      int offset = ml.offsets[i];
      int size = ml.offsets[i + 1] - ml.offsets[i];
      if (offset >= 0 && size > 0)
      {
        int b_id = start_buffer_offset + i;
        clear_buffer(b_id, 0.0f);
        for (int j = 0; j < vertex_count; j++)
          memcpy(buffers[b_id] + size*j, part_data + ml.f_per_vert * j + offset, sizeof(float)*size);
        set_array_to_ctx_internal(buffer_names[b_id], b_id, size * vertex_count);
      }
    }
    show_errors();
    start_buffer_offset += 3;
  }
}

void MitsubaInterface::camera_to_ctx(const CameraSettings &camera, std::string camera_name)
{
  PyObject *func, *args, *p_name, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *p9, *p10, *p11, *p12, *func_ret;

  func = PyObject_GetAttrString(pModule, (char *)"set_camera");

  p_name = PyUnicode_FromString(camera_name.c_str());

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

  args = PyTuple_Pack(14, mitsubaContext, p_name, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12);
  func_ret = PyObject_CallObject(func, args);
  show_errors();

  DEL(func);
  DEL(args);
  DEL(p_name);
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

void MitsubaInterface::render_model_to_file(const dgen::DFModel &model, const std::string &image_dir,
                                            const CameraSettings &camera, const std::vector<float> &scene_params)
{  
  model_to_ctx(model);
  camera_to_ctx(camera, "camera");

  int cameras_buf_id = get_camera_buffer_id();
  assert(scene_params.size() > 0);
  memcpy(buffers[cameras_buf_id], scene_params.data(), sizeof(float)*scene_params.size());
  set_array_to_ctx_internal(buffer_names[cameras_buf_id], cameras_buf_id, scene_params.size());
  show_errors();

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

float MitsubaInterface::render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, const std::vector<float> &scene_params,
                                           double *timers)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  model_to_ctx(model);
  for (int i=0;i<cameras.size();i++)
    camera_to_ctx(cameras[i], "camera_"+std::to_string(i));

  int cameras_buf_id = get_camera_buffer_id();
  assert(scene_params.size() > 0);
  memcpy(buffers[cameras_buf_id], scene_params.data(), sizeof(float)*scene_params.size());
  set_array_to_ctx_internal(buffer_names[cameras_buf_id], cameras_buf_id, scene_params.size());
  show_errors();

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  float loss = render_and_compare_internal();
  std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
  //get derivatives by vertex positions and scene parameters
  get_array_from_ctx_internal(buffer_names[0] + "_grad", 0);
  get_array_from_ctx_internal(buffer_names[cameras_buf_id] + "_grad", cameras_buf_id);
  std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
  if (timers)
  {
    timers[2] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    timers[3] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    timers[4] += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
  }
  return loss;
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
  if (data_floats > 4*model_max_size)
  {
    logerr("Python array %s contains %d float, while buffer size is %d. Some data will be ignored", name.c_str(), data_floats, 4*model_max_size);
  }
  char *data = PyBytes_AsString(params_bytes);
  memcpy(buffers[buffer_id], data, MIN(sz, 4*model_max_size*sizeof(float)));
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

void MitsubaInterface::render_multicam_demo(RenderSettings render_settings, ModelInfo model_info, 
                                            const dgen::DFModel &model, const std::string &image_dir,
                                            const std::vector<float> &scene_params, const CameraSettings &camera,
                                            int rotations_x, int rotations_y)
{
  auto &ml = model_info.layout;
  auto &rs = render_settings;

  //find center of model to adjust cameras
  glm::vec3 center = glm::vec3(0,0,0);
  for (int i=0;i<model.first.size();i+=ml.f_per_vert)
  {
    center += glm::vec3(model.first[i + ml.pos], model.first[i + ml.pos + 1], model.first[i + ml.pos + 2]);
  }
  center = center/(float)(model.first.size()/ml.f_per_vert);
  float dist = glm::length(center - camera.origin);

  Texture composite_tex = engine::textureManager->create_texture(rotations_x*rs.image_w, 2*rotations_y*rs.image_h);
  PostFx copy("copy.fs");
  int fbo = create_framebuffer();

  auto render_to_tex = [&](glm::ivec2 offset)
  {
    for (int i=0;i<rotations_x;i++)
    {
      for (int j=0;j<rotations_y;j++)
      {
        float phi = (2*PI*i)/rotations_x;
        float psi = (2*PI*j)/(rotations_y + 1);//no need for top view

        glm::vec3 pos = dist*glm::vec3(cos(psi)*cos(phi), sin(psi), cos(psi)*sin(phi));
        glm::vec3 right = glm::cross(pos, glm::vec3(0,1,0));
        glm::vec3 up = glm::normalize(glm::cross(right, pos));
        CameraSettings cur_cam = camera;
        cur_cam.target = center;
        cur_cam.origin = center + pos;
        cur_cam.up = up;

        render_model_to_file(model, image_dir, cur_cam, scene_params);
        Texture raw = engine::textureManager->load_unnamed_tex(image_dir, 1);
        engine::view->next_frame();
        
        int prev_FBO = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, composite_tex.texture, 0);
        glViewport(offset.x + i*rs.image_w, offset.y + j*rs.image_h, rs.image_w, rs.image_h);
        copy.use();
        copy.get_shader().texture("tex", raw);
        copy.render();
        
        glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
      }
    }
  };

  //Render with constant texture
  RenderSettings setting_monochrome(rs.image_w, rs.image_h, rs.samples_per_pixel, rs.mitsubaVar, RenderStyle::MONOCHROME_DEMO);
  init_scene_and_settings(setting_monochrome, model_info);
  render_to_tex(glm::ivec2(0,0));

  RenderSettings setting_textured(rs.image_w, rs.image_h, rs.samples_per_pixel, rs.mitsubaVar, RenderStyle::TEXTURED_DEMO);
  init_scene_and_settings(setting_textured, model_info);
  render_to_tex(glm::ivec2(0,rotations_y*rs.image_h));

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  engine::textureManager->save_png_directly(composite_tex, image_dir);

  delete_framebuffer(fbo);
}

std::vector<std::string> MitsubaInterface::get_all_available_materials()
{
  return 
  {
    "very smooth porcelain",
    "smooth porcelain",
    "porcelain",
    "ceramics",
    "rough ceramics",
    "glass",
    "imperfect glass",
    "frosted glass",
    "lambert"
  };
}

std::string MitsubaInterface::get_default_material()
{
  return "ceramics";
}

CameraSettings MitsubaInterface::get_camera_from_scene_params(const std::vector<float> &scene_params)
{
  float fov_rad = scene_params.back();

  CameraSettings camera;
  float h1 = 1.5;
  camera.fov_rad = fov_rad;
  float h2 = h1 * tan((3.14159265f / 3) / 2) / tan(camera.fov_rad / 2);
  camera.origin = glm::vec3(0, 0.5, h2);
  camera.target = glm::vec3(0, 0.5, 0);
  camera.up = glm::vec3(0, 1, 0);
  return camera;
}