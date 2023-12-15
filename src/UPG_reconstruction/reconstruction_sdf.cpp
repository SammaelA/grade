#include "upg.h"
#include "sdf_node.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "common_utils/bbox.h"
#include "sdf_rendering.h"
#include "graphics_utils/render_point_cloud.h"

namespace upg
{

  struct PointCloudReference
  {
    std::vector<glm::vec3> points;
    std::vector<glm::vec3> outside_points;
    bool is_synthetic = false;
    UPGStructure structure;//can be set manually to make reconstruction simplier
    UPGParametersRaw parameters;//empty if not synthetic reference
  };

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
      sdf_to_point_cloud(sdf, points, &(reference.points), &(reference.outside_points));

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

    // estimate sum(f(0), ..., f(max_index)) by sampling only part of these values
    // function expect all f(i) >= 0 and increases number of samples if average f(i) is small
    // returns partial sum and number of samples
    static std::pair<int, double> sum_with_adaptive_batching(std::function<double(int)> get_value, 
                                                             int max_index, int base_batch_size,
                                                             int max_batch_size,
                                                             double sensitivity)
    {
      int samples = 0;
      double partial_sum = 0.0;
      while (samples < max_batch_size && partial_sum <= sensitivity)
      {
        for (int b=0;b<base_batch_size;b++)
        {
          int index = rand()%max_index;
          partial_sum += get_value(index);
        }
        samples += base_batch_size;
        sensitivity /= 10;
        if (sensitivity == 0)
          break;
      }

      return {samples, partial_sum};
    }

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(reference.points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);
      assert(gen_params.size() == out_grad.size());

      std::vector<float> cur_grad;
      std::vector<float> dpos_dparams = {0,0,0};
      cur_grad.reserve(gen_params.size());
      std::vector<double> out_grad_d(gen_params.size(), 0);

      //main step - minimize SDF values on surface
      auto p1 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.points[index];
        cur_grad.clear();
        double d = sdf.get_distance(p, &cur_grad, &dpos_dparams);
        for (int i=0;i<gen_params.size();i++)
          out_grad_d[i] += 2*d*cur_grad[i];
        return d*d;
      }, reference.points.size(), 256, 5000, 10);

      //regularization - penalty for outside points with sdf < 0
      auto p2 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.outside_points[index];
        cur_grad.clear();
        double d = sdf.get_distance(p, &cur_grad, &dpos_dparams);
        if (d < 0)
        {
          for (int i=0;i<gen_params.size();i++)
            out_grad_d[i] += 2*d*cur_grad[i];
          return d*d;
        }
        else
          return 0;
      }, reference.outside_points.size(), 128, 5000, 0.1);

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
      
      for (int i=0;i<gen_params.size();i++)
        out_grad[i] = out_grad_d[i]/(p1.first+p2.first);
      
      return p1.second/p1.first + p2.second/p2.first;
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(reference.points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);

      //main step - minimize SDF values on surface
      auto p1 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.points[index];
        double d = sdf.get_distance(p);
        return d*d;
      }, reference.points.size(), 256, 5000, 10);

      //regularization - penalty for outside points with sdf < 0
      auto p2 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.outside_points[index];
        double d = sdf.get_distance(p);
        if (d < 0)
          return d*d;
        else
          return 0;
      }, reference.outside_points.size(), 128, 5000, 0.1);

      return p1.second/p1.first + p2.second/p2.first;
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
      else if (optimizer_name == "particle_swarm")
        optimizer = get_optimizer_particle_swarm(&opt_func, *step_blk, start_params.structure);
      else if (optimizer_name == "CC")
        optimizer = get_optimizer_CC(&opt_func, *step_blk, start_params.structure);
      else if (optimizer_name == "DE")
        optimizer = get_optimizer_differentiable_evolution(&opt_func, *step_blk, start_params.structure);
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

      PointCloudRenderer renderer;
      Texture tp = renderer.render(reference.points, camera.get_viewProj(), 1024, 1024, {1,0,0}, 0.2);
      engine::textureManager->save_png(tp, "result_points");
    }

    return opt_res;
  }
}