#include "upg.h"
#include "sdf_node.h"
#include "common_utils/distribution.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "common_utils/bbox.h"
#include "preprocessing.h"
#include "graphics_utils/image_metrics.h"

namespace upg
{
  inline glm::vec3 EyeRayDirNormalized(float x/*in [0,1]*/, float y/*in [0,1]*/, glm::mat4 projInv)
  {
    glm::vec4 pos = glm::vec4(2.0f*x - 1.0f, -2.0f*y + 1.0f, 0.0f, 1.0f);
    pos = projInv * pos;
    pos /= pos.w;
    return glm::normalize(glm::vec3(pos));
  }

  inline glm::vec3 transformRay(glm::vec3 ray, glm::mat4 viewInv)
  {
    glm::vec3 p1 = glm::vec3(viewInv*glm::vec4(0,0,0,1));
    glm::vec3 p2 = glm::vec3(viewInv*glm::vec4(ray.x,ray.y,ray.z,1));
    return glm::normalize(p2-p1);
  }

  Texture render_sdf(const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, bool lambert = true)
  {
    glm::mat4 projInv = glm::inverse(camera.get_proj());
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    //set light somewhere to the side 
    glm::vec3 light_dir = normalize(camera.origin + glm::vec3(camera.origin.z, camera.origin.y, camera.origin.x) - camera.target);
    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    #pragma omp parallel for
    for (int yi=0;yi<image_h;yi++)
    {
      std::vector<float> cur_grad;
      std::vector<float> ddist_dpos = {0,0,0};
      for (int xi=0;xi<image_w;xi++)
      {
        glm::vec3 color = {0,0,0};
        for (int yp=0;yp<spp_a;yp++)
        {
          for (int xp=0;xp<spp_a;xp++)
          {
            float y = (float)(yi*spp_a+yp)/(image_h*spp_a);
            float x = (float)(xi*spp_a+xp)/(image_w*spp_a);
            glm::vec3 p0 = camera.origin;
            glm::vec3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            
            //sphere tracing
            int iter = 0;
            float d = sdf.get_distance(p0);
            while (iter < 1000 && d > 1e-6 && d < 1e6)
            {
              p0 += d*dir;
              d = sdf.get_distance(p0);
              iter++;
            }
            if (d <= 1e-6)
            {
              if (lambert)
              {
                cur_grad.clear();
                sdf.get_distance(p0, &cur_grad, &ddist_dpos);
                glm::vec3 n = glm::normalize(glm::vec3(ddist_dpos[0], ddist_dpos[1], ddist_dpos[2]));
                color += glm::vec3(1,1,1) * MAX(0.1f, dot(n, light_dir));
              }
              else
                color += glm::vec3(1,1,1);
            }
          }
        }
        if (!lambert)
          color = {1,1,1};
        data[4*(yi*image_w+xi)+0] = 255*(color.x/SQR(spp_a));
        data[4*(yi*image_w+xi)+1] = 255*(color.y/SQR(spp_a));
        data[4*(yi*image_w+xi)+2] = 255*(color.z/SQR(spp_a));
        data[4*(yi*image_w+xi)+3] = 255;
      }
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  struct PointCloudReference
  {
    std::vector<glm::vec3> points;
    std::vector<glm::vec3> outside_points;
    bool is_synthetic = false;
    UPGStructure structure;//can be set manually to make reconstruction simplier
    UPGParametersRaw parameters;//empty if not synthetic reference
  };

  AABB get_point_cloud_bbox(const std::vector<glm::vec3> &points);
  void sdf_to_point_cloud(const ProceduralSdf &sdf, int points, PointCloudReference &cloud)
  {
    float r = 10000;

    cloud.points.reserve(points);
    for (int i=0;i<points;i++)
    {
      float phi = urand(0, 2*PI);
      float psi = urand(-PI/2, PI/2);
      glm::vec3 p0 = r*glm::vec3{cos(psi)*cos(phi), sin(psi), cos(psi)*sin(phi)};
      glm::vec3 dir = glm::normalize(-p0); //tracing rays from random points to (0,0,0)

      int iter = 0;
      float d = sdf.get_distance(p0);
      while (iter < 1000 && d > 1e-6 && d < 1e6)
      {
        p0 += d*dir;
        d = sdf.get_distance(p0);
        iter++;
      }
      if (d <= 1e-6)
        cloud.points.push_back(p0);
    }

    AABB bbox = get_point_cloud_bbox(cloud.points);
    AABB inflated_bbox = AABB(bbox.min_pos - glm::vec3(0.01,0.01,0.01), bbox.max_pos + glm::vec3(0.01,0.01,0.01));
    
    cloud.outside_points.reserve(points);
    while (cloud.outside_points.size() < points)
    {
      glm::vec3 p = glm::vec3(urand(inflated_bbox.min_pos.x, inflated_bbox.max_pos.x),
                              urand(inflated_bbox.min_pos.y, inflated_bbox.max_pos.y),
                              urand(inflated_bbox.min_pos.z, inflated_bbox.max_pos.z));
      if (!bbox.contains(p))
        cloud.outside_points.push_back(p);
    }
  }

  AABB get_point_cloud_bbox(const std::vector<glm::vec3> &points)
  {
    glm::vec3 minv(1e9,1e9,1e9);
    glm::vec3 maxv(-1e9,-1e9,-1e9);
    for (auto &p : points)
    {
      minv = glm::min(minv, p);
      maxv = glm::max(maxv, p);
    }

    return AABB(minv, maxv);
  }

  float get_sdf_image_based_quality(ProceduralSdf reference_sdf, ProceduralSdf sdf)
  {
    int image_size = 512;
    int points = 1000;
    PointCloudReference cloud;
    sdf_to_point_cloud(reference_sdf, points, cloud);
    AABB bbox = get_point_cloud_bbox(cloud.points);
    float d = 1.5*length(bbox.max_pos - bbox.min_pos);
    CameraSettings cam;
    cam.target = 0.5f*(bbox.min_pos + bbox.max_pos);
    cam.origin = cam.target + glm::vec3(0,0,-d);

    auto cameras = get_cameras_uniform_sphere(cam, 64, d);

    float diff = 0.0;
    float mse = 0.0;
    for (auto &cam : cameras)
    {
      Texture t1 = render_sdf(reference_sdf, cam, image_size, image_size, 4);
      Texture t2 = render_sdf(sdf, cam, image_size, image_size, 4);
      mse += ImageMetric::get(t1, t2);
    }
    return -10*log10(MAX(1e-9f,mse/cameras.size()));
  }

  PointCloudReference get_point_cloud_reference(const Block &input_blk)
  {
    PointCloudReference reference;
    Block *synthetic_reference = input_blk.get_block("synthetic_reference"); //reference is parameters for our own generator
    Block *model_reference = input_blk.get_block("model_reference"); //reference is an .obj (or other) file with 3D model
    assert(!(synthetic_reference && model_reference));
    if (synthetic_reference)
    {
      reference.is_synthetic = true;
      synthetic_reference->get_arr("structure", reference.structure.s);
      synthetic_reference->get_arr("params", reference.parameters.p);
      SdfGenInstance gen(reference.structure);
      ProceduralSdf sdf = gen.generate(reference.parameters.p);

      int points = synthetic_reference->get_int("points_count", 10000);
      sdf_to_point_cloud(sdf, points, reference);

      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 16);
      engine::textureManager->save_png(t, "reference_sdf");
    }
    else if (model_reference)
    {
      //TODO
    }

    return reference;
  };

  class SdfRenderAndCompare : public UPGOptimizableFunction
  {
  public:
    SdfRenderAndCompare(const PointCloudReference &_reference, const Block &optimization_blk):
    reference(_reference)
    {
      
    }

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(reference.points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);

      std::vector<float> cur_grad;
      std::vector<float> dpos_dparams = {0,0,0};
      cur_grad.reserve(gen_params.size());
      for (int i=0;i<gen_params.size();i++)
        out_grad[i] = 0;

      double full_d = 0.0;

      //main step - minimize SDF values on surface
      for (const glm::vec3 &p : reference.points)
      {
        cur_grad.clear();
        float d = sdf.get_distance(p, &cur_grad, &dpos_dparams);
        full_d += SQR(d);
        for (int i=0;i<gen_params.size();i++)
          out_grad[i] += 2*d*cur_grad[i] / reference.points.size();
      }

      for (const glm::vec3 &p : reference.outside_points)
      {
        cur_grad.clear();
        float d = sdf.get_distance(p, &cur_grad, &dpos_dparams) - 0.01;
        if (d < 0)
        {
          full_d += SQR(d);
          for (int i=0;i<gen_params.size();i++)
            out_grad[i] += 2*d*cur_grad[i] / reference.points.size();
        }
      }
      //regularization - penalty for outside points with sdf < 0
      /*
      debug("params [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",params.differentiable[i]);
      debug("]\n");
      debug("grad [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",out_grad[i]);
      debug("]\n");
      */
      return full_d/reference.points.size();
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(reference.points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);
      double d = 0.0;
      for (const glm::vec3 &p : reference.points)
        d += SQR(sdf.get_distance(p));
      for (const glm::vec3 &p : reference.outside_points)
      {
        float od = sdf.get_distance(p) - 0.01;
        if (od < 0)
          d += SQR(od);
      }

      return d/reference.points.size();
    }
    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) override
    {
      return ((SdfGenInstance*)gen)->desc;
    }
    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const override
    {
      return std::make_shared<SdfGenInstance>(structure);
    }
  private:
    const PointCloudReference &reference;
  };

  std::vector<UPGReconstructionResult> reconstruct_sdf(const Block &blk)
  {
    //load settings from given blk
    Block *input_blk = blk.get_block("input");
    Block *gen_blk   = blk.get_block("generator");
    Block *opt_blk   = blk.get_block("optimization");
    Block *res_blk   = blk.get_block("results");
    if (!input_blk || !gen_blk || !opt_blk || !res_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    //get ReconstructionReference - all info about the object that we want to reconstruct
    PointCloudReference reference = get_point_cloud_reference(*input_blk);
    SdfGenInstance::set_scene_bbox(get_point_cloud_bbox(reference.points));

    //get start parameters for optimization. They are required for Adam and other local optimizers
    //and have to be set manually
    UPGReconstructionResult start_params;
    Block *start_params_blk = opt_blk->get_block("start");
    if (start_params_blk)
    {
      start_params_blk->get_arr("params", start_params.parameters.p);
      start_params_blk->get_arr("structure", start_params.structure.s);
    }

    //perform optimization. There might be one or several steps of it, I expect the first step 
    //to be some sort of Genetic Algorithm and others - Adam optimizers for fine-tuning the params
    int step_n = 0;
    std::vector<UPGReconstructionResult> opt_res = {start_params};
    while (opt_blk->get_block("step_"+std::to_string(step_n)))
    {
      Block *step_blk = opt_blk->get_block("step_"+std::to_string(step_n));

      SdfRenderAndCompare opt_func(reference, *step_blk);
      std::shared_ptr<UPGOptimizer> optimizer;
      std::string optimizer_name = step_blk->get_string("optimizer_name", "adam");
      if (optimizer_name == "adam")
        optimizer = get_optimizer_adam(&opt_func, *step_blk, opt_res[0]);
      else if (optimizer_name == "memetic")
        optimizer = get_optimizer_memetic(&opt_func, *step_blk, start_params.structure);
      else if (optimizer_name == "CHC")
        optimizer = get_optimizer_CHC(&opt_func, *step_blk, start_params.structure);
      optimizer->optimize();
      opt_res = optimizer->get_best_results();
      step_n++;
    }

    for (auto &result : opt_res)
    {
      SdfGenInstance gen(result.structure);
      ProceduralSdf sdf = gen.generate(result.parameters.p);

      result.quality_ir = result.loss_optimizer;

      if (reference.is_synthetic && res_blk->get_bool("check_model_quality"))
      {
        SdfGenInstance reference_gen(reference.structure);
        ProceduralSdf reference_sdf = reference_gen.generate(reference.parameters.p);
        result.quality_synt = get_sdf_image_based_quality(reference_sdf, sdf);
      }

      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 16);
      engine::textureManager->save_png(t, "result_sdf");
    }

    return opt_res;
  }
}