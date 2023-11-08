#include "upg.h"
#include "preprocessing.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"
#include "tinyEngine/engine.h"
#include "custom_diff_render/custom_diff_render.h"
#include "simple_render_and_compare.h"
#include "graphics_utils/modeling.h"
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
    UPGOptimizer(const Block &settings, ReconstructionReference &reference)
    {
      render_w = settings.get_int("render_w", 128);
      render_h = settings.get_int("render_h", 128);
      cameras_pd = get_cameras_parameter_description(reference.images);
      diff_render.reset(get_halfgpu_custom_diff_render());
      simple_render.reset(new NonDiffRender());
      
      //get diff_render settings
      IDiffRender::Settings diff_render_settings;
      diff_render_settings.image_w = render_w;
      diff_render_settings.image_h = render_h;

      std::vector<Texture> references;
      for (int i=0; i<reference.images.size(); i++)
      {
        reference.images[i].resized_mask = resize_mask(reference.images[i].mask, render_w, render_h, true);
        references.push_back(reference.images[i].resized_mask);
      }

      if (diff_render)
      diff_render->init_optimization(references, diff_render_settings, 
                                     settings.get_bool("save_intermediate_images", false));
      if (simple_render)
      simple_render->init_optimization(references, diff_render_settings, 
                                     settings.get_bool("save_intermediate_images", false));
    }
  protected:
    //all parameters that can be changed by optimizer structured
    //in a convenient (or optimizer) way
    struct Params
    {
      std::vector<float> differentiable;
    };

    void opt_params_to_gen_params_and_camera(const Params &params, const ParametersDescription &pd, 
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
      opt_params_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenMesh mesh = gen.generate(full_gen_params);
      UniversalGenJacobian dPos_dP = gen.generate_jacobian(full_gen_params);
      float res = diff_render->render_and_compare_silhouette(mesh.pos, cameras);
      std::span<const float> dLoss_dPos(diff_render->get_vertex_grad(), mesh.pos.size());
      grad_jac_mult(dPos_dP, dLoss_dPos, out_grad);

      return res;
    }

    int render_w, render_h;
    std::unique_ptr<IDiffRender> diff_render;
    std::unique_ptr<IDiffRender> simple_render;
    ParametersDescription cameras_pd;
  };

  class UPGOptimizerAdam : public UPGOptimizer
  {
  public:
    UPGOptimizerAdam(const Block &settings, ReconstructionReference &reference, const UPGReconstructionResult &start_params) :
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
          debug("Adam iter %3d  val = %.4f best_val = %.4f\n", iter, val, best_result);
      }
      if (verbose)
        debug("Adam final res val = %.4f best_val = %.4f\n", best_result, best_result);

      UPGReconstructionResult res;
      res.structure = gen_structure;
      {
        std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
        std::vector<CameraSettings> cameras;
        opt_params_to_gen_params_and_camera(best_params, pd, full_gen_params, cameras);
        res.parameters.p = full_gen_params;
        res.loss_optimizer = best_result;
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
    //load settings from given blk
    Block *input_blk = blk.get_block("input");
    Block *gen_blk = blk.get_block("generator");
    Block *opt_blk = blk.get_block("optimization");
    Block *res_blk = blk.get_block("results");
    if (!input_blk || !gen_blk || !opt_blk || !res_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    //get ReconstructionReference - all info about the object that we want to reconstruct
    ReconstructionReference reference = get_reference(*input_blk);

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
      std::unique_ptr<UPGOptimizer> optimizer(new UPGOptimizerAdam(*step_blk, reference, opt_res[0]));
      opt_res = optimizer->optimize();
      step_n++;
    }
    
    //TODO: compare results with reference and calculate reconstruction quality
    //Also calculate another quality metric if we have synthetic reference
    for (auto &result : opt_res)
    {
      ComplexModel reconstructed_model;
      std::string save_directory = "saves/" + res_blk->get_string("save_folder") + "/";
      if (res_blk->get_string("save_folder") != "")
        prepare_directory("saves/" + res_blk->get_string("save_folder"));
      if (!create_model(result.structure, result.parameters, reconstructed_model))
        logerr("failed to create model from reconstructed structure and parameters!");
      if (res_blk->get_bool("check_image_quality"))
        result.quality_ir = get_image_based_quality(reference, reconstructed_model);
      if (reference.model.is_valid() && res_blk->get_bool("check_model_quality"))
        result.quality_synt = get_model_based_quality(reference, reconstructed_model);
    
      if (res_blk->get_bool("save_turntable") && res_blk->get_block("save_turntable_hydra_settings"))
        render_model_turntable(*(res_blk->get_block("save_turntable_hydra_settings")), reconstructed_model);
      if (reference.model.is_valid() && res_blk->get_bool("save_reference_turntable") && 
          res_blk->get_block("save_reference_turntable_hydra_settings"))
        render_model_turntable(*(res_blk->get_block("save_reference_turntable_hydra_settings")), reference.model);

      if (res_blk->get_bool("save_model"))
      {
        assert(reconstructed_model.is_valid());
        assert(reconstructed_model.models.size() == 1); //TODO: save complex models with different materials too
        model_loader::save_model_to_obj(reconstructed_model.models[0], save_directory + "reconstructed_model.obj");

        if (reference.model.is_valid())
        {
          assert(reference.model.models.size() == 1); //TODO: save complex models with different materials too
          model_loader::save_model_to_obj(reference.model.models[0], save_directory + "reference_model.obj");          
        }
      } 
    }

    //print results of the reconstruction
    if (opt_blk->get_bool("verbose"))
    {
      debug("UPG Reconstruction finished\n");
      for (int i=0;i<opt_res.size();i++)
      {
        debug("=============================\n");
        debug("Optimization result %d\n",i);
        debug("Optimizer's loss    :%7.5f\n", opt_res[i].loss_optimizer);
        debug("Image-based quality :%7.5f\n", opt_res[i].quality_ir);
        if (reference.model.is_valid())
        debug("Model-based quality :%7.5f\n", opt_res[i].quality_synt);
        else
        debug("Model-based quality : ----- \n");
      }
      debug("=============================\n");
    }
    
    return opt_res;
  }
}