#include "upg.h"
#include "sdf_node.h"
#include "common_utils/distribution.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"

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

  Texture render_sdf(const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp)
  {
    glm::mat4 projInv = glm::inverse(camera.get_proj());
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    for (int yi=0;yi<image_h;yi++)
    {
      for (int xi=0;xi<image_w;xi++)
      {
        float c = 0;
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
            while (iter < 100 && d > 1e-6 && d < 1e6)
            {
              p0 += d*dir;
              d = sdf.get_distance(p0);
              iter++;
            }
            if (d <= 1e-6)
              c++;
          }
        }
        unsigned char color = 255*(c/SQR(spp_a));
        data[4*(yi*image_w+xi)+0] = color;
        data[4*(yi*image_w+xi)+1] = color;
        data[4*(yi*image_w+xi)+2] = color;
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
  };

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
      while (iter < 100 && d > 1e-6 && d < 1e6)
      {
        p0 += d*dir;
        d = sdf.get_distance(p0);
        iter++;
      }
      if (d <= 1e-6)
        cloud.points.push_back(p0);
    }
  }

  PointCloudReference get_point_cloud_reference(const Block &input_blk)
  {
    PointCloudReference reference;
    Block *synthetic_reference = input_blk.get_block("synthetic_reference"); //reference is parameters for our own generator
    Block *model_reference = input_blk.get_block("model_reference"); //reference is an .obj (or other) file with 3D model
    assert(!(synthetic_reference && model_reference));
    if (synthetic_reference)
    {
      UPGStructure structure;
      UPGParametersRaw params;
      synthetic_reference->get_arr("structure", structure.s);
      synthetic_reference->get_arr("params", params.p);
      SdfGenInstance gen(structure);
      ProceduralSdf sdf = gen.generate(params.p);

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
      for (const glm::vec3 &p : reference.points)
      {
        cur_grad.clear();
        float d = SQR(sdf.get_distance(p, &cur_grad, &dpos_dparams));
        full_d += SQR(d);
        for (int i=0;i<gen_params.size();i++)
          out_grad[i] += 2*d*cur_grad[i] / reference.points.size();
      }
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
      optimizer->optimize();
      opt_res = optimizer->get_best_results();
      step_n++;
    }

    return opt_res;
  }
}