#include "preprocessing.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/postfx.h"
#include "graphics_utils/simple_model_utils.h"
#include "simple_render_and_compare.h"
#include "core/scene.h"
#include "hydra/hydra_scene_exporter.h"
#include "graphics_utils/modeling.h"
#include <random>
#include <glm/glm.hpp>

namespace upg
{
  Texture render_silhouette(const ComplexModel &model, const CameraSettings &camera, int tex_w, int tex_h)
  {
    Shader render_silhouette({"render_silhouette.vs", "render_silhouette.fs"}, {"in_Position", "in_Normal", "in_Tex"});
    Texture tex = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);

    int fbo = create_framebuffer();
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex.texture, 0);
    glViewport(0, 0, tex_w, tex_h);

    glm::mat4 y_swap = glm::mat4(glm::vec4(1,0,0,0), glm::vec4(0,-1,0,0), glm::vec4(0,0,1,0),glm::vec4(0,0,0,1));
    //render model silhouette
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glm::mat4 projection = glm::perspective(camera.fov_rad, 1.0f, camera.z_near, camera.z_far);
    glm::mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);
    glm::mat4 viewProj = y_swap * projection * view;
    render_silhouette.use();
    render_silhouette.uniform("viewProj", viewProj);
    for (auto &m : model.models)
      m->render();

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
    delete_framebuffer(fbo);

    return tex;
  }

  Texture get_mask(Texture &image)
  {
    // TODO: image to mask
    logerr("preprocessing: get_mask is not implemented");
    return {};
  }

  ReferenceView get_reference_view(const Block &view_blk)
  {
    ReferenceView rv;
    rv.camera = visualizer::load_camera_settings(view_blk);
    rv.fixed_camera = view_blk.get_bool("camera.fixed", true);

    std::string mask_dir = view_blk.get_string("mask", "");
    std::string image_dir = view_blk.get_string("image", "");

    if (mask_dir == "")
    {
      if (image_dir == "")
        logerr("preprocessing: each view block must have mask or image path");
      Texture image = engine::textureManager->load_unnamed_tex(image_dir, 1);
      rv.mask = get_mask(image);
    }
    else
    {
      rv.mask = engine::textureManager->load_unnamed_tex(mask_dir, 1);
    }

    return rv;
  }

  ReferenceView get_reference_view_from_model(const ComplexModel &model, const Block &model_reference_blk, const Block &view_blk)
  {
    ReferenceView rv;
    rv.camera = visualizer::load_camera_settings(view_blk);
    rv.fixed_camera = view_blk.get_bool("camera.fixed", true);
    rv.mask = render_silhouette(model, rv.camera, 
                                model_reference_blk.get_int("reference_image_w", 1024),
                                model_reference_blk.get_int("reference_image_h", 1024));

    return rv;
  }

  ReconstructionReference get_reference(const Block &input_blk)
  {
    ReconstructionReference reference;
    Block *synthetic_reference = input_blk.get_block("synthetic_reference"); //reference is parameters for our own generator
    Block *model_reference = input_blk.get_block("model_reference"); //reference is an .obj (or other) file with 3D model
    assert(!(synthetic_reference && model_reference));
    if (synthetic_reference)
    {
      create_model_from_block(*synthetic_reference, reference.model);
      reference.model.update();
    }
    else if (model_reference)
    {
      //load obj
      std::string obj_filename = model_reference->get_string("obj_filename");
      std::string tex_filename = model_reference->get_string("tex_filename");
      assert(obj_filename != "");
      Texture t;
      if (tex_filename == "") 
        t = engine::textureManager->get("porcelain");
      else
        t = engine::textureManager->load_unnamed_tex(tex_filename,1);
      Model *m = model_loader::load_model_from_obj_directly(obj_filename);
      if (model_reference->get_bool("normalize_model"))
        model_loader::normalize_model(m);
      reference.model.models = {m};
      reference.model.materials = {Material(t)};
      reference.model.update();
    }

    for (int i = 0; i < input_blk.size(); i++)
    {
      Block *view_blk = input_blk.get_block(std::string("view_")+std::to_string(i));
      if (view_blk)
      {
        if (synthetic_reference)
          reference.images.push_back(get_reference_view_from_model(reference.model, *synthetic_reference, *view_blk));
        else if (model_reference)
          reference.images.push_back(get_reference_view_from_model(reference.model, *model_reference, *view_blk));
        else
          reference.images.push_back(get_reference_view(*view_blk));

        //engine::textureManager->save_png(reference.images.back().mask, "reference_"+std::to_string(i));
      }
    }

    return reference;
  }

  Texture resize_mask(Texture mask, int tex_w, int tex_h, bool to_binary_mask)
  {
    Texture resized_mask = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);
    PostFx copy("copy.fs");
    PostFx binarize("binarize.fs");

    int fbo = create_framebuffer();
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, resized_mask.texture, 0);
    glViewport(0, 0, tex_w, tex_h);

    glClear(GL_COLOR_BUFFER_BIT);
    copy.use();
    copy.get_shader().texture("tex", mask);
    copy.render();
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    if (to_binary_mask)
    {
      Texture binary_mask = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, binary_mask.texture, 0);
      glViewport(0, 0, tex_w, tex_h);

      glClear(GL_COLOR_BUFFER_BIT);
      binarize.use();
      binarize.get_shader().texture("tex", resized_mask);
      binarize.render();
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
      
      resized_mask = binary_mask;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
    delete_framebuffer(fbo);

    return resized_mask;
  }

  float get_image_based_quality(const ReconstructionReference &reference, const ComplexModel &reconstructed_model)
  {
    NonDiffRender render;
    IDiffRender::Settings render_settings;
    render_settings.image_w = reference.images[0].mask.get_W();
    render_settings.image_h = reference.images[0].mask.get_H();
    std::vector<CameraSettings> cameras;
    std::vector<Texture> images;
    for (auto &im : reference.images)
    {
      assert(render_settings.image_w == im.mask.get_W());
      assert(render_settings.image_h == im.mask.get_H());
      cameras.push_back(im.camera);
      images.push_back(im.mask);
    }
    render.init_optimization(images, render_settings);
    
    std::vector<float> positions;
    for (auto *m : reconstructed_model.models)
      positions.insert(positions.end(), m->positions.begin(), m->positions.end());
    
    float mse = render.render_and_compare_silhouette(positions, cameras);
    return -10*log10(MAX(1e-9f,mse));
  }

  std::vector<CameraSettings> get_cameras_uniform_sphere(CameraSettings orig_camera, int cams_n, float distance)
  {
    std::vector<CameraSettings> cameras;
    std::mt19937 gen(0);
    std::uniform_real_distribution<float> d_ur = std::uniform_real_distribution<float>(0,1);

    for (int i=0;i<cams_n;i++)
    {
      float phi = 2*PI*d_ur(gen);
      float psi = (PI/2)*(1 - sqrtf(d_ur(gen)));
      if (d_ur(gen) > 0.5)
        psi = - psi;
      
      cameras.push_back(orig_camera);
      glm::vec3 view_dir = -glm::vec3(cos(psi)*sin(phi), sin(psi), cos(psi)*cos(phi));
      glm::vec3 tangent = glm::normalize(glm::cross(view_dir, orig_camera.up));
      glm::vec3 new_up = glm::normalize(glm::cross(view_dir, tangent));
      cameras.back().up = new_up;
      cameras.back().origin = cameras.back().target -distance*view_dir;
    }

    return cameras;
  }

  float get_model_based_quality(const ReconstructionReference &reference, const ComplexModel &reconstructed_model)
  {
    if (!reference.model.is_valid())
      return 0;
    
    NonDiffRender render;
    IDiffRender::Settings render_settings;
    render_settings.image_w = reference.images[0].mask.get_W();
    render_settings.image_h = reference.images[0].mask.get_H();
    float distance = length(reference.images[0].camera.origin - reference.images[0].camera.target);
    std::vector<CameraSettings> cameras = get_cameras_uniform_sphere(reference.images[0].camera, 64, distance);
    std::vector<Texture> images;
    for (auto &c : cameras)
      images.push_back(render_silhouette(reference.model, c, render_settings.image_w, render_settings.image_h));
    render.init_optimization(images, render_settings);
    
    std::vector<float> positions;
    for (auto *m : reconstructed_model.models)
      positions.insert(positions.end(), m->positions.begin(), m->positions.end());
    
    float mse = render.render_and_compare_silhouette(reconstructed_model.models[0]->positions, cameras);
    return -10*log10(MAX(1e-9f,mse));
  }

  void render_model_turntable(const Block &hydra_settings, const ComplexModel &model)
  {
    Scene scene;
    scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
    scene.heightmap->fill_const(0);
    scene.instanced_models.emplace_back();
    scene.instanced_models[0].instances = {glm::mat4(1.0f)};
    scene.instanced_models[0].model = model;
    Block export_settings;
    hydra::prepare_hydra_export_settings_block(hydra_settings, export_settings);
    hydra::export_scene("model_turntable", scene, export_settings);
  }
}