#include "preprocessing.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/shader.h"
#include "graphics_utils/simple_model_utils.h"
#include <glm/glm.hpp>

namespace upg
{
  Texture render_silhouette(const ComplexModel &model, const CameraSettings &camera, int tex_w, int tex_h)
  {
    Shader render_silhouette({"render_silhouette.vs", "render_silhouette.fs"}, {"in_Position", "in_Normal", "in_Tex"});
    Texture tex = engine::textureManager->create_texture(tex_w, tex_h, GL_RGB8);

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
    engine::textureManager->save_png(tex, "preprocessing_test.png");
    return tex;
  }

  Texture get_mask(Texture &image)
  {
    // TODO: image to mask
    logerr("preprocessing: get_mask is not implemented");
    return {};
  }

  ReferenceView preprocess_get_reference_view(const Block &view_blk)
  {
    ReferenceView rv;
    rv.camera = visualizer::load_camera_settings(view_blk);
    rv.fixed_camera = view_blk.get_bool("camera.fixed", true);

    Block *synt_reference = view_blk.get_block("synthetic_reference");
    if (synt_reference)
    {
      // our reference is a set of procedural generator's parameters
      // we should create model from these parameters and render it
      // with given camera (it contains in view_blk too)
      ComplexModel model;
      create_model_from_block(*synt_reference, model);
      model.update();
      rv.mask = render_silhouette(model, rv.camera, 
                                  synt_reference->get_int("tex_w", 1024),
                                  synt_reference->get_int("tex_h", 1024));
    }
    else
    {
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
    }

    return rv;
  }
}