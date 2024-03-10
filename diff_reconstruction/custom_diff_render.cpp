#include "custom_diff_render.h"
#include "common_utils/utility.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/modeling.h"
#include "common_utils/LiteMath_ext.h"
#include <chrono>

#ifdef USE_CUSTOM_DIFF_RENDER
#include "../modules/custom_diff_render/src/drender.h"
#include "LiteMath/Image2d.h"
class CustomDiffRender : public IDiffRender
{
public:
  CustomDiffRender()
  {
    
  }
  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images = false) override
  {
    if (render_settings.renderStyle != MitsubaInterface::SILHOUETTE)
      logerr("CustomDiffRender supports only SILHOUETTE render style");

    diff_render::DiffRenderSettings settings{diff_render::SHADING_MODEL::SILHOUETTE, render_settings.samples_per_pixel};
    _dr = diff_render::MakeDifferentialRenderer(settings);

    for (auto &ref : reference_image_dir)
      reference_images.push_back(LiteImage::LoadImage<diff_render::float3>(ref.c_str()));
    
    save_images = save_intermediate_images;
    if (save_images)
    {
      for (auto &ref : reference_image_dir)
        out_images.push_back(diff_render::Img(reference_images[0].width(), reference_images[0].height()));
    }
  }

  virtual float render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params) override
  {
    auto t1 = std::chrono::steady_clock::now();
    //convert model to diff_render::Scene
    diff_render::Scene scene;
    diff_render::TriangleMesh mesh;
    int cnt = model.first.size()/FLOAT_PER_VERTEX;
    mesh.vertices.resize(cnt);
    mesh.tc.resize(cnt, {0,0});
    mesh.normals.resize(cnt, {1,0,0});//we don't need normals and tc for silhouette rendering
    mesh.indices.resize(cnt);
    for (int i=0;i<cnt;i++)
    {
      mesh.vertices[i] = {model.first[i*FLOAT_PER_VERTEX], model.first[i*FLOAT_PER_VERTEX+1], model.first[i*FLOAT_PER_VERTEX+2]};
      mesh.indices[i] = i;
    }

    //TODO fix me. The same strange transform applied in Mitsuba!
    float tr_z = scene_params[2] * tan(0.5*0.25) / tan(0.5*cameras[0].fov_rad);

    diff_render::TransformR tr({scene_params[0], scene_params[1], tr_z}, {scene_params[3], scene_params[4], scene_params[5]}, 1);
    scene.add_mesh(mesh, {tr});

    //create DScene
    scene_grad = diff_render::DScene(scene, diff_render::SHADING_MODEL::SILHOUETTE, {0});

    //convert cameras to diff_render::CamInfo
    std::vector<diff_render::CamInfo> dr_cameras;
    int cam_w = reference_images[0].width();
    int cam_h = reference_images[0].height();
    for (auto &cam : cameras)
      dr_cameras.push_back(diff_render::CamInfo({cam.origin.x, cam.origin.y, cam.origin.z},
                                                {cam.target.x, cam.target.y, cam.target.z},
                                                {cam.up.x, cam.up.y, cam.up.z}, 
                                                cam_w, cam_h, cam.fov_rad, cam.z_near, cam.z_far));

    //render to get gradients
    assert(dr_cameras.size() == reference_images.size());
    auto t2 = std::chrono::steady_clock::now();
    float loss = _dr->d_render_and_compare(scene, dr_cameras.data(), reference_images.data(), dr_cameras.size(), cam_w*cam_h, scene_grad,
                                           save_images ? out_images.data() : nullptr);
    auto t3 = std::chrono::steady_clock::now();
    //store scene_parameters_grad separately, as we can calculate grad only for first scene params = transform+rotate;
    scene_parameters_grad = std::vector<float>(scene_params.size(), 0);
    for (int i=0;i<6;i++)
      scene_parameters_grad[i] = scene_grad.get_dmeshes()[0].restricted_transform(0)[i];
  
    if (save_images)
    {
      for (int i=0; i<out_images.size(); i++)
      {
        std::string img_name = "output/img_"+std::to_string(iteration)+"_cam_"+std::to_string(i)+".png";
        LiteImage::SaveImage(img_name.c_str(), out_images[i]);
      }
    }

    auto t4 = std::chrono::steady_clock::now();
    int t_all = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t1).count();
    int t_real = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    logerr("Custom diff. render took = %d (%d + %d) ms", t_all, t_real, t_all - t_real);

    iteration++;
    static float average_ms = 0.0;
    average_ms += t_all;
    logerr("average val %f", average_ms/(iteration));
    return 10*log10(MAX(1e-9f,loss));
  }
  
  virtual const float *get_vertex_grad() const override
  {
    assert(sizeof(float) == sizeof(diff_render::GradReal));
    return scene_grad.get_dmeshes()[0].pos(0);
  }

  virtual const float *get_scene_params_grad() const override
  {
    return scene_parameters_grad.data();
  }

private:
  std::shared_ptr<diff_render::IDiffRender> _dr;
  std::vector<diff_render::Img> reference_images;
  std::vector<diff_render::Img> out_images;
  mutable diff_render::DScene scene_grad;
  std::vector<float> scene_parameters_grad;
  bool save_images = false;
  int iteration = 0;

};

class HalfGPUCustomDiffRender : public IDiffRender
{
public:
  HalfGPUCustomDiffRender():
  render_silhouette({"render_silhouette.vs", "render_silhouette.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
  diff_loss("diff_loss.fs"),
  diff_loss_sum({"diff_loss_sum.comp"},{})
  {
    fbo = create_framebuffer();

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    
    results_buf = create_buffer();

    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
  }

  virtual ~HalfGPUCustomDiffRender()
  {
    delete_framebuffer(fbo);
    delete_buffer(results_buf);
  }

  virtual void init_optimization(const std::vector<std::string> &reference_image_dir, 
                                MitsubaInterface::RenderSettings render_settings, 
                                bool save_intermediate_images = false) override
  {
    if (render_settings.renderStyle != MitsubaInterface::SILHOUETTE)
      logerr("HalfGPUCustomDiffRender supports only SILHOUETTE render style");

    diff_render::DiffRenderSettings settings{diff_render::SHADING_MODEL::SILHOUETTE, render_settings.samples_per_pixel};
    _dr = diff_render::MakeDifferentialRenderer(settings);

    for (auto &ref : reference_image_dir)
    {
      reference_images.push_back(LiteImage::LoadImage<diff_render::float3>(ref.c_str()));
      adjoints.push_back(diff_render::Img(reference_images[0].width(), reference_images[0].height()));
    }
    
    save_images = save_intermediate_images;
    if (save_images)
    {
      for (auto &ref : reference_image_dir)
        out_images.push_back(diff_render::Img(reference_images[0].width(), reference_images[0].height()));
    }

 
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  assert(reference_image_dir.size() >= 1);
  if (render_settings.renderStyle != MitsubaInterface::SILHOUETTE)
    logerr("DiffRenderGPU supports only silhouette render style");

  tex_w = render_settings.image_w;
  tex_h = render_settings.image_h;
  save_images = save_intermediate_images;
  int edge_samples = tex_w*tex_h;

  for (auto &dir : reference_image_dir)
  {
    reference_textures.push_back(engine::textureManager->load_unnamed_tex(dir, 1));
    assert(reference_textures.back().get_W() == tex_w && reference_textures.back().get_H() == tex_h);
  }

  tex1 = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);
  tex2 = engine::textureManager->create_texture(tex_w, tex_h, GL_RGB16F);

  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
 
  }

  virtual float render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras, 
                                  const std::vector<float> &scene_params) override
  {
    auto t0 = std::chrono::steady_clock::now();

    assert(cameras.size() == reference_textures.size());
    Model m;
    visualizer::simple_mesh_to_model_332(model.first, &m);
    m.update();

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex1.texture, 0);
    glViewport(0, 0, tex_w, tex_h);

    //TODO fix me. The same strange transform applied in Mitsuba!
    float tr_z = scene_params[2] * tan(0.5*0.25) / tan(0.5*cameras[0].fov_rad);
    glm::mat4 transform = LiteMath::translate(glm::mat4(1.0f),glm::vec3(scene_params[0], scene_params[1], tr_z))*
                          LiteMath::rotate(glm::mat4(1.0f), scene_params[5], glm::vec3(0,0,1))*
                          LiteMath::rotate(glm::mat4(1.0f), scene_params[4], glm::vec3(0,1,0))*
                          LiteMath::rotate(glm::mat4(1.0f), scene_params[3], glm::vec3(1,0,0));
    float total_res = 0;

    for (int i=0;i<cameras.size();i++)
    {
      //render model silhouette
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glm::mat4 viewProj = cameras[i].get_viewProj();
      render_silhouette.use();
      render_silhouette.uniform("viewProj", viewProj);
      m.render();

      //create diff_loss texture (0.5*(1 + img - target_img))
      glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex2.texture, 0);
      diff_loss.use();
      diff_loss.get_shader().uniform("tex_size", glm::vec2(tex_w, tex_h));
      diff_loss.get_shader().texture("tex", tex1);
      diff_loss.get_shader().texture("tex_reference", reference_textures[i]);
      diff_loss.render();

      //save rendered image if needed
      if (save_images)
      {
        glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
        std::string dir = std::string("DiffRenderGPU_")+std::to_string(iteration)+std::string("_cam_")+std::to_string(i);
        engine::textureManager->save_png(tex1, dir);
      }

      #define NUM_TILES 16

      //calculate MSE per tile
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
      diff_loss_sum.use();
      diff_loss_sum.uniform("tex_size", glm::vec2(tex_w, tex_h));
      diff_loss_sum.texture("tex_diff", tex2);
      glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
      glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*NUM_TILES*NUM_TILES, nullptr, GL_STREAM_READ);
      glDispatchCompute(1, 1, 1);

      //accumulate per-tile MSE to one value
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
      float res = 0;
      float results[NUM_TILES*NUM_TILES];
      glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
      GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
      memcpy(results,ptr,sizeof(float)*NUM_TILES*NUM_TILES);
      for (int j=0;j<NUM_TILES*NUM_TILES;j++)
        res += results[j];
      total_res += res;
      glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);  

      //read back texture
      {
        Texture &t = tex2;
        glBindTexture(t.type, t.texture);

        glGetTexImage(t.type,
                      0,
                      GL_RGB,
                      GL_FLOAT,
                      adjoints[i].data());
        glBindTexture(t.type, 0);
      }
    }

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

    auto t1 = std::chrono::steady_clock::now();
    //convert model to diff_render::Scene
    diff_render::Scene scene;
    diff_render::TriangleMesh mesh;
    int cnt = model.first.size()/FLOAT_PER_VERTEX;
    mesh.vertices.resize(cnt);
    mesh.tc.resize(cnt, {0,0});
    mesh.normals.resize(cnt, {1,0,0});//we don't need normals and tc for silhouette rendering
    mesh.indices.resize(cnt);
    for (int i=0;i<cnt;i++)
    {
      mesh.vertices[i] = {model.first[i*FLOAT_PER_VERTEX], model.first[i*FLOAT_PER_VERTEX+1], model.first[i*FLOAT_PER_VERTEX+2]};
      mesh.indices[i] = i;
    }

    diff_render::TransformR tr({scene_params[0], scene_params[1], tr_z}, {scene_params[3], scene_params[4], scene_params[5]}, 1);
    scene.add_mesh(mesh, {tr});

    //create DScene
    scene_grad = diff_render::DScene(scene, diff_render::SHADING_MODEL::SILHOUETTE, {0});

    //convert cameras to diff_render::CamInfo
    std::vector<diff_render::CamInfo> dr_cameras;
    int cam_w = reference_images[0].width();
    int cam_h = reference_images[0].height();
    for (auto &cam : cameras)
      dr_cameras.push_back(diff_render::CamInfo({cam.origin.x, cam.origin.y, cam.origin.z},
                                                {cam.target.x, cam.target.y, cam.target.z},
                                                {cam.up.x, cam.up.y, cam.up.z}, 
                                                cam_w, cam_h, cam.fov_rad, cam.z_near, cam.z_far));

    //render to get gradients
    assert(dr_cameras.size() == reference_images.size());
    auto t2 = std::chrono::steady_clock::now();    
    _dr->commit(scene);
    scene_grad.clear();
    _dr->d_render(scene, dr_cameras.data(), adjoints.data(), dr_cameras.size(), cam_w*cam_h, scene_grad);
    auto t3 = std::chrono::steady_clock::now();

    //store scene_parameters_grad separately, as we can calculate grad only for first scene params = transform+rotate;
    scene_parameters_grad = std::vector<float>(scene_params.size(), 0);
    for (int i=0;i<6;i++)
      scene_parameters_grad[i] = scene_grad.get_dmeshes()[0].restricted_transform(0)[i];
  
    if (save_images)
    {
      for (int i=0; i<out_images.size(); i++)
      {
        std::string img_name = "output/img_"+std::to_string(iteration)+"_cam_"+std::to_string(i)+".png";
        LiteImage::SaveImage(img_name.c_str(), out_images[i]);
      }
    }

    auto t4 = std::chrono::steady_clock::now();
    int t_all = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t0).count();
    int t_real = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    int t_render = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    //logerr("Custom diff. render took = %.1f (%.1f + %.1f + %.1f) ms", 1e-3*t_all, 1e-3*t_render, 1e-3*t_real, 1e-3*(t_all - t_real - t_render));

    iteration++;
    static float average_ms = 0.0;
    average_ms += 1e-3*t_all;
    //logerr("average val %f", average_ms/(iteration));

    return 10*log10(MAX(1e-9f,total_res/cameras.size()));
  }
  
  virtual const float *get_vertex_grad() const override
  {
    assert(sizeof(float) == sizeof(diff_render::GradReal));
    return scene_grad.get_dmeshes()[0].pos(0);
  }

  virtual const float *get_scene_params_grad() const override
  {
    return scene_parameters_grad.data();
  }

private:
  std::shared_ptr<diff_render::IDiffRender> _dr;
  std::vector<diff_render::Img> reference_images;
  std::vector<diff_render::Img> adjoints;
  std::vector<diff_render::Img> out_images;
  mutable diff_render::DScene scene_grad;
  std::vector<float> scene_parameters_grad;
  bool save_images = false;
  int iteration = 0;

  GLuint fbo, results_buf;
  Shader render_silhouette, diff_loss_sum;
  PostFx diff_loss;
  int tex_w = 0, tex_h = 0;
  std::vector<Texture> reference_textures;
  Texture tex1;
  Texture tex2;

};

IDiffRender *create_custom_diff_render()
{
  return new HalfGPUCustomDiffRender();
}
#else
int custom_diff_render_main(int argc, char *argv[])
{
  logerr("custom differentiable renderer is not implemented!");
  return 0;
}
IDiffRender *create_custom_diff_render()
{
  logerr("custom differentiable renderer is not implemented!");
  return nullptr;
}
#endif