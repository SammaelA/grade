#include "reconstruction.h"
#include "preprocessing.h"
#include "reconstruction_impl.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"
#include "tinyEngine/engine.h"
#include "custom_diff_render/custom_diff_render.h"
#include <memory>
#include <unistd.h>

namespace upg
{
  std::vector<ParametersDescription::Param> get_camera_params(const CameraSettings &cam)
  {
    std::vector<ParametersDescription::Param> params;

    params.push_back({cam.origin.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.x"});
    params.push_back({cam.origin.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.y"});
    params.push_back({cam.origin.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.z"});

    params.push_back({cam.target.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.x"});
    params.push_back({cam.target.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.y"});
    params.push_back({cam.target.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.z"});

    params.push_back({cam.up.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.x"});
    params.push_back({cam.up.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.y"});
    params.push_back({cam.up.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.z"});

    params.push_back({cam.z_near, 0.001f, 1.0f, ParameterType::MUTABLE_FLOAT, "z_near"});
    params.push_back({cam.z_far, 10.0f, 1000.0f, ParameterType::MUTABLE_FLOAT, "z_far"});
    params.push_back({cam.fov_rad, 0.1f, 1.5f, ParameterType::MUTABLE_FLOAT, "fov_rad"});

    return params;
  }

  unsigned get_camera_block_id(int camera_n)
  {
    return (1 << 16) + (unsigned)camera_n;
  }

  ParametersDescription get_cameras_parameter_description(const std::vector<ReferenceView> &reference_views)
  {
    ParametersDescription pd;
    for (int i = 0; i < reference_views.size(); i++)
    {
      auto pv = get_camera_params(reference_views[i].camera);
      if (reference_views[i].fixed_camera)
      {
        for (auto &p : pv)
          p.type = ParameterType::CONST;
      }
      pd.add_parameters(get_camera_block_id(i), "camera_" + std::to_string(i), pv);
    }
    return pd;
  }

  std::vector<ReferenceView> get_reference(const Block &input_blk)
  {
    std::vector<ReferenceView> reference;
    for (int i = 0; i < input_blk.size(); i++)
    {
      Block *view_blk = input_blk.get_block(i);
      if (!view_blk)
      {
        logerr("UPG Reconstruction: views block should contain only blocks with data for each view");
        continue;
      }
      reference.push_back(preprocess_get_reference_view(*view_blk));
    }

    return reference;
  }

  Texture resize_mask(Texture mask, int tex_w, int tex_h)
  {
    assert(mask.get_W() == tex_w && mask.get_H() == tex_h);
    return mask;
  }
  CameraSettings camera_from_params(const std::vector<float> &p)
  {
    CameraSettings cam;
    cam.origin = glm::vec3(p[0], p[1], p[2]);
    cam.target = glm::vec3(p[3], p[4], p[5]);
    cam.up     = glm::vec3(p[6], p[7], p[8]);
    cam.z_near = p[9];
    cam.z_far = p[10];
    cam.fov_rad = p[11];
    return cam;
  }

  
  void grad_jac_mult(const UniversalGenJacobian &dPos_dP, std::span<const float> dLoss_dPos, 
                     std::span<float> out_dLoss_dP)
  {
    assert(out_dLoss_dP.size() == dPos_dP.y_n);
    assert(dLoss_dPos.size() == dPos_dP.x_n);
    for (int i=0;i<dPos_dP.y_n;i++)
    { out_dLoss_dP[i] = 0;
      for (int j=0;j<dPos_dP.x_n;j++)
        out_dLoss_dP[i] += dPos_dP.jacobian[i*dPos_dP.x_n + j]*dLoss_dPos[j]; 
    }
  }

  class UPGOptimizer
  {
  public:
    virtual ~UPGOptimizer() = default;
    virtual std::vector<UPGReconstructionResult> optimize() = 0;
    UPGOptimizer(const Block &settings, std::vector<ReferenceView> &reference)
    {
      render_w = settings.get_int("render_w", 128);
      render_h = settings.get_int("render_h", 128);
      cameras_pd = get_cameras_parameter_description(reference);
      diff_render.reset(get_halfgpu_custom_diff_render());
      
      //get diff_render settings
      IDiffRender::Settings diff_render_settings;
      diff_render_settings.image_w = render_w;
      diff_render_settings.image_h = render_h;

      std::vector<std::string> reference_paths;
      for (int i=0; i<reference.size(); i++)
      {
        reference[i].resized_mask = resize_mask(reference[i].mask, render_w, render_h);
        std::string reference_path = "saves/reference_"+std::to_string(i)+".png";
        reference_paths.push_back(reference_path);
        engine::textureManager->save_png_directly(reference[i].mask, reference_path);
      }
      sleep(1); //to be sure that png save is finished

      diff_render->init_optimization(reference_paths, diff_render_settings, 
                                     settings.get_bool("save_intermediate_images", false));
    }
  protected:
    //all parameters that can be changed by optimizer structured
    //in a convenient (or optimizer) way
    struct Params
    {
      std::vector<float> differentiable;
    };

    void opt_paras_to_gen_params_and_camera(const Params &params, const ParametersDescription &pd, 
                                            /*out*/ std::vector<float> &full_gen_params,
                                            /*out*/ std::vector<CameraSettings> &cameras)
    {
      //iterate all parameters' groups from description
      //map orders them by block_id, so generator's params are first, and cameras' after it
      int diff_i = 0;
      full_gen_params.reserve(pd.get_total_params_count());
      for (const auto &p : pd.get_block_params())
      {
        std::vector<float> camera_params;
        bool is_camera_block = p.first >= get_camera_block_id(0);
        std::vector<float> &p_v = is_camera_block ? camera_params : full_gen_params;
        for (auto &par : p.second.p)
        {
          if (par.type == ParameterType::CONST)
            p_v.push_back(par.value);
          else if (par.type == ParameterType::DIFFERENTIABLE)
          {
            float v = CLAMP(params.differentiable[diff_i], par.min_val, par.max_val);
            p_v.push_back(v);
            diff_i++;
          }
          else
          {
            // TODO: other types
          }
        }
        if (is_camera_block)
          cameras.push_back(camera_from_params(camera_params));
      }
    }

    //calculate function that we optimize and it's gradient (put into given span)
    //requires already created UniversalGenInstance
    //size of out_grad - is a number of differentiable parameters in params
    float f_grad_f(UniversalGenInstance &gen, const ParametersDescription &pd,
                   const Params &params, std::span<float> out_grad)
    {
      std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
      std::vector<CameraSettings> cameras;
      opt_paras_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenMesh mesh = gen.generate(full_gen_params);
      UniversalGenJacobian dPos_dP = gen.generate_jacobian(full_gen_params);
      float res = diff_render->render_and_compare_silhouette(mesh.pos, cameras);
      std::span<const float> dLoss_dPos(diff_render->get_vertex_grad(), mesh.pos.size());
      grad_jac_mult(dPos_dP, dLoss_dPos, out_grad);

      return res;
    }

    int render_w, render_h;
    std::unique_ptr<IDiffRender> diff_render;
    ParametersDescription cameras_pd;
  };

  class UPGOptimizerAdam : public UPGOptimizer
  {
  public:
    UPGOptimizerAdam(const Block &settings, std::vector<ReferenceView> &reference, const UPGReconstructionResult &start_params) :
    UPGOptimizer(settings, reference),
    gen(start_params.structure)
    {
      iterations = settings.get_int("iterations");
      alpha = settings.get_double("learning_rate", 0.01);
      beta_1 = settings.get_double("beta_1", 0.9);
      beta_2 = settings.get_double("beta_2", 0.999);
      eps = settings.get_double("eps", 1e-8);
      verbose = settings.get_bool("verbose") || settings.get_int("verbose") > 0;

      X_n = start_params.parameters.p.size();
      X.differentiable = start_params.parameters.p;
      gen_structure = start_params.structure;

      pd.add(cameras_pd);
      pd.add(gen.desc);
      pd.print_info();
    }
    virtual std::vector<UPGReconstructionResult> optimize() override
    {
      std::vector<float> V = std::vector<float>(X_n, 0); 
      std::vector<float> S = std::vector<float>(X_n, 0);
      UPGOptimizer::Params best_params = X;
      std::vector<float> x_grad = std::vector<float>(X_n, 0); 
      float best_result = 1e9;

      for (int iter=0; iter<iterations; iter++)
      {
        float val = f_grad_f(gen, pd, X, x_grad);
        if (val < best_result)
        {
          best_params = X;
          best_result = val;
        }
        for (int i=0;i<X_n;i++)
        {
          float g = x_grad[i];
          V[i] = beta_1 * V[i] + (1-beta_1)*g;
          float Vh = V[i] / (1 - pow(beta_1, iter+1)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
          float Sh = S[i] / (1 - pow(beta_2, iter+1)); 
          X.differentiable[i] -= alpha*Vh/(sqrt(Sh) + eps);
        }
        if ((iter % 5 == 0) && verbose)
          debug("Adam iter %d val = %.4f best_val = %.4f\n", iter, val, best_result);
      }

      UPGReconstructionResult res;
      res.structure = gen_structure;
      {
        std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
        std::vector<CameraSettings> cameras;
        opt_paras_to_gen_params_and_camera(best_params, pd, full_gen_params, cameras);
        res.parameters.p = full_gen_params;
        res.quality = -10*log10(MAX(1e-9f,best_result)); //MSE to PSNR
      }
      return {res};
    }

  protected:
    int iterations;
    float alpha, beta_1, beta_2, eps;
    bool verbose;
    int X_n;
    UPGOptimizer::Params X;
    UniversalGenInstance gen;
    ParametersDescription pd;
    UPGStructure gen_structure;
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk)
  {
    //we store settings in separate block for each stage of the algorithm
    Block *input_blk = blk.get_block("input");
    Block *gen_blk = blk.get_block("generator");
    Block *opt_blk = blk.get_block("optimization");
    if (!input_blk || !gen_blk || !opt_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    std::vector<ReferenceView> reference = get_reference(*input_blk);
    UPGReconstructionResult start_params;
    start_params.parameters.p = {0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81};

    std::unique_ptr<UPGOptimizer> optimizer(new UPGOptimizerAdam(*opt_blk, reference, start_params));
    auto opt_res = optimizer->optimize();
    return opt_res;
  }
}