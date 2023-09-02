#include "custom_diff_render.h"
#include "common_utils/utility.h"
#include <chrono>

#ifdef USE_CUSTOM_DIFF_RENDER
#include "drender.h"
#include "Image2d.h"
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
IDiffRender *create_custom_diff_render()
{
  return new CustomDiffRender();
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