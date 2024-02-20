#include "upg.h"
#include "sdf_node.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "common_utils/bbox.h"
#include "sdf_rendering.h"
#include "graphics_utils/render_point_cloud.h"
#include "common_utils/distribution.h"
#include "sdf_reconstruction_common.h"
#include <chrono>

namespace upg
{
  class PointCloudSdfLoss : public UPGOptimizableFunction
  {
  public:
    PointCloudSdfLoss(const PointCloudReference &_reference, const Block &optimization_blk):
    reference(_reference)
    {
      
    }

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(reference.points.size() > 0);
      //if we don't have constant parameters opt_params and gen_params are exact same vector
      //and there is no need to copy them. Except for the clamping min/max values
      std::vector<float> explicit_gen_params;
      if (true)
        explicit_gen_params = opt_params_to_gen_params(params, pd);
      const std::vector<float> &gen_params = true ? explicit_gen_params : params.data;
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);
      assert(gen_params.size() == out_grad.size());

      if (positions.size() < 3*max_batch_size)
        positions = std::vector<float>(3*max_batch_size,0);
      if (distances.size() < max_batch_size)
        distances = std::vector<float>(max_batch_size,0);
      if (dparams.size() < gen_params.size()*max_batch_size)
        dparams = std::vector<float>(gen_params.size()*max_batch_size,0);
      if (dpositions.size() < 3*max_batch_size)
        dpositions = std::vector<float>(3*max_batch_size,0);

      unsigned batch_size = 256;
      unsigned param_count = gen_params.size();
      std::vector<double> out_grad_d(param_count, 0);
      double loss = 0;
      double reg_loss = 0;

      //main step - minimize SDF values on surface
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.points.size();
        positions[3*i+0] = reference.points[index].x;
        positions[3*i+1] = reference.points[index].y;
        positions[3*i+2] = reference.points[index].z;
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), dparams.data(), dpositions.data());

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = glm::sign(distances[i])*std::min(0.03f, abs(distances[i]));
        for (int j=0;j<param_count;j++)
          out_grad_d[j] += 2*d*dparams[i*param_count + j];
        loss += d*d;
      }

      //regularization - penalty for outside points with sdf < 0
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.outside_points.size();
        positions[3*i+0] = reference.outside_points[index].x;
        positions[3*i+1] = reference.outside_points[index].y;
        positions[3*i+2] = reference.outside_points[index].z;
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), dparams.data(), dpositions.data());

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = distances[i];
        if (d < 0)
        {
          for (int j=0;j<param_count;j++)
            out_grad_d[j] += (d > 0 ? 2 : 4)*d*dparams[i*param_count + j];
          reg_loss += d*d;
        }
      }

      for (int i=0;i<gen_params.size();i++)
        out_grad[i] = out_grad_d[i]/(2*batch_size);
      /* debug("params1 [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",params[i]);
      debug("]\n");
      debug("grad1 [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",out_grad[i]);
      debug("]\n");  */   

      return loss/batch_size + reg_loss/batch_size;
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(reference.points.size() > 0);
      //if we don't have constant parameters opt_params and gen_params are exact same vector
      //and there is no need to copy them
      std::vector<float> explicit_gen_params;
      if (pd.has_constants())
        explicit_gen_params = opt_params_to_gen_params(params, pd);
      const std::vector<float> &gen_params = pd.has_constants() ? explicit_gen_params : params.data;
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);

      if (positions.size() < 3*max_batch_size)
        positions = std::vector<float>(3*max_batch_size,0);
      if (distances.size() < max_batch_size)
        distances = std::vector<float>(max_batch_size,0);

      unsigned batch_size = 256;
      unsigned param_count = gen_params.size();
      double loss = 0;
      double reg_loss = 0;

      //main step - minimize SDF values on surface
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.points.size();
        positions[3*i+0] = reference.points[index].x;
        positions[3*i+1] = reference.points[index].y;
        positions[3*i+2] = reference.points[index].z;
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), nullptr, nullptr);

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = glm::sign(distances[i])*std::min(0.03f, abs(distances[i]));
        loss += d*d;
      }

      //regularization - penalty for outside points with sdf < 0
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.outside_points.size();
        positions[3*i+0] = reference.outside_points[index].x;
        positions[3*i+1] = reference.outside_points[index].y;
        positions[3*i+2] = reference.outside_points[index].z;
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), nullptr, nullptr);

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = distances[i];
        if (d < 0)
          reg_loss += d*d;
      } 
      //logerr("loss %f %f\n",(float)loss/batch_size, (float)reg_loss/batch_size);
      return loss/batch_size + reg_loss/batch_size;
    }
    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) const override
    {
      return ((ProceduralSdf*)gen)->desc;
    }
    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const override
    {
      return std::make_shared<ProceduralSdf>(structure);
    }

    virtual float estimate_positioning_quality(const UPGStructure &structure,
                                               const UPGPart &part, std::span<const float> parameters,
                                               float border_sigma,
                                               float inner_point_penalty) const override
  {
    UPGStructure part_structure;
    part_structure.s = std::vector<uint16_t>(structure.s.begin() + part.s_range.first, structure.s.begin() + part.s_range.second);
    std::span<const float> part_parameters(parameters.data() + part.p_range.first, parameters.data() + part.p_range.second);

    ProceduralSdf sdf(part_structure);
    sdf.set_parameters(part_parameters);

    double quality = 0, q1 = 0, q2 = 0;
    float denom = normal_pdf(0,0,border_sigma);
    for (auto &p : reference.points)
    {
      float dist = sdf.get_distance(p);
      float border_q = normal_pdf(dist, 0, border_sigma)/denom; //in range (0,1], where 1 is perfect border point, 0 is not border at all
      float inner_q = (dist<0)*(1-border_q);
      quality += border_q - inner_point_penalty*inner_q;
      q1 += border_q;
      q2 += inner_q;
    }
    quality /= reference.points.size();

    //logerr("part [%d %d][%d %d]", (int)part.s_range.first, (int)part.s_range.second, part.p_range.first, part.p_range.second);
    //logerr("quality %f (%f %f)", (float)quality, (float)q1, (float)q2);

    return quality;
  }
  private:
    const PointCloudReference &reference;
    unsigned max_batch_size = 1024;
    std::vector<float> positions, distances, dparams, dpositions;
  };

  std::vector<UPGReconstructionResult> simple_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                  const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    PointCloudSdfLoss opt_func(reference, *step_blk);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer(step_blk->get_string("optimizer_name", "adam"), &opt_func, step_blk,
                                                            prev_step_res[0], prev_step_res[0].structure);
    optimizer->optimize();
    return optimizer->get_best_results();
  }

  std::vector<UPGReconstructionResult> reconstruction_graph_based(Block *step_blk, const std::vector<glm::vec3> &points,
                                                                  const std::vector<float> &distances);
  std::vector<UPGReconstructionResult> constructive_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                        const std::vector<UPGReconstructionResult> &prev_step_res);

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
    ProceduralSdf::set_scene_bbox(get_point_cloud_bbox(reference.points));

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
    auto t1 = std::chrono::steady_clock::now();
    while (opt_blk->get_block("step_"+std::to_string(step_n)))
    {
      Block *step_blk = opt_blk->get_block("step_"+std::to_string(step_n));
      if (step_blk->get_bool("constructive_reconstruction"))
        opt_res = constructive_reconstruction_step(step_blk, reference, opt_res);
      else if (step_blk->get_bool("graph_based"))
        opt_res = reconstruction_graph_based(step_blk, reference.d_points, reference.d_distances);
      else
        opt_res = simple_reconstruction_step(step_blk, reference, opt_res);
      step_n++;
    }
    auto t2 = std::chrono::steady_clock::now();
    float t_ms = 0.001f*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logerr("reconstruction took %.3f s",t_ms/1000);

    for (auto &result : opt_res)
    {
      ProceduralSdf sdf(result.structure);
      sdf.set_parameters(result.parameters.p);

      result.quality_ir = result.loss_optimizer;

      if (reference.is_synthetic && res_blk->get_bool("check_model_quality"))
      {
        ProceduralSdf reference_sdf(reference.structure);
        reference_sdf.set_parameters(reference.parameters.p);
        result.quality_ir = get_sdf_similarity_MSE(reference_sdf, sdf);
      }
      if (reference.is_synthetic && res_blk->get_bool("check_image_quality"))
      {
        ProceduralSdf reference_sdf(reference.structure);
        reference_sdf.set_parameters(reference.parameters.p);
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
      //engine::textureManager->save_png(tp, "result_points");
    }

    return opt_res;
  }
}