#include "upg.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"
#include "tinyEngine/engine.h"
#include "custom_diff_render/custom_diff_render.h"
#include "simple_render_and_compare.h"
#include "graphics_utils/modeling.h"
#include "optimization.h"
#include <memory>
#include <unistd.h>
#include <algorithm>

namespace upg
{
  class MeshRenderAndCompare : public UPGOptimizableFunction
  {
  public:
    MeshRenderAndCompare(const ReconstructionReference &reference, const Block &optimization_blk):
    diff_render(get_halfgpu_custom_diff_render()),
    simple_render(new NonDiffRender())
    {

      int render_w = optimization_blk.get_int("render_w", 128);
      int render_h = optimization_blk.get_int("render_h", 128);
      cameras_pd = get_cameras_parameter_description(reference.images);
      
      //get diff_render settings
      IDiffRender::Settings diff_render_settings;
      diff_render_settings.image_w = render_w;
      diff_render_settings.image_h = render_h;

      std::vector<Texture> references;
      for (int i=0; i<reference.images.size(); i++)
        references.push_back(resize_mask(reference.images[i].mask, render_w, render_h, true));

      if (diff_render)
      diff_render->init_optimization(references, diff_render_settings, 
                                     optimization_blk.get_bool("save_intermediate_images", false));
      if (simple_render)
      simple_render->init_optimization(references, diff_render_settings, 
                                       optimization_blk.get_bool("save_intermediate_images", false));
    }

    virtual float f_grad_f(UniversalGenInstance &gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
      std::vector<CameraSettings> cameras;
      opt_params_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenJacobian dPos_dP;
      UniversalGenMesh mesh = gen.generate(full_gen_params, &dPos_dP);
      float res = diff_render->render_and_compare_silhouette(mesh.pos, cameras);
      std::span<const float> dLoss_dPos(diff_render->get_vertex_grad(), mesh.pos.size());
      grad_jac_mult(dPos_dP, dLoss_dPos, out_grad);

      return res;
    }

    virtual float f_no_grad(UniversalGenInstance &gen, const ParametersDescription &pd, const OptParams &params) override
    {
      std::vector<float> full_gen_params; //all parameters, including non-differentiable and consts, for generation. No cameras here
      std::vector<CameraSettings> cameras;
      opt_params_to_gen_params_and_camera(params, pd, full_gen_params, cameras);
      UniversalGenMesh mesh = gen.generate(full_gen_params, nullptr);
      
      return simple_render->render_and_compare_silhouette(mesh.pos, cameras);
    }
    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance &gen) override
    {
      ParametersDescription pd;
      pd.add(gen.desc);
      pd.add(cameras_pd);
      return pd;
    }

  private:
    static void grad_jac_mult(const UniversalGenJacobian &dPos_dP, std::span<const float> dLoss_dPos, 
                              std::span<float> out_dLoss_dP)
    {
      assert(out_dLoss_dP.size() == dPos_dP.get_yn());
      assert(dLoss_dPos.size() == dPos_dP.get_xn());
      for (int i=0;i<dPos_dP.get_yn();i++)
      { out_dLoss_dP[i] = 0;
        for (int j=0;j<dPos_dP.get_xn();j++)
          out_dLoss_dP[i] += dPos_dP.at(i,j)*dLoss_dPos[j]; 
      }
    }
    static std::vector<ParametersDescription::Param> get_camera_params(const CameraSettings &cam)
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

    static unsigned get_camera_block_id(int camera_n)
    {
      return get_scene_params_block_offset() + (unsigned)camera_n;
    }

    static ParametersDescription get_cameras_parameter_description(const std::vector<ReferenceView> &reference_views)
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

    static CameraSettings camera_from_params(const std::span<float> &p)
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
    void opt_params_to_gen_params_and_camera(const OptParams &params, const ParametersDescription &pd, 
                                             /*out*/ std::vector<float> &full_gen_params,
                                             /*out*/ std::vector<CameraSettings> &cameras)
    {
      std::vector<float> camera_params;
      full_gen_params = opt_params_to_gen_params(params, pd, &camera_params);
      constexpr int num_camera_params = 12;
      assert(cameras.empty());
      assert(camera_params.size() > 0 && (camera_params.size() % num_camera_params == 0)); 

      for (int i=0;i<camera_params.size();i+=num_camera_params)
        cameras.push_back(camera_from_params(std::span<float>(camera_params.data() + i, num_camera_params)));
    };

    std::unique_ptr<IDiffRender> diff_render;
    std::unique_ptr<IDiffRender> simple_render;
    ParametersDescription cameras_pd;
  
  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk)
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

      MeshRenderAndCompare opt_func(reference, *step_blk);
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