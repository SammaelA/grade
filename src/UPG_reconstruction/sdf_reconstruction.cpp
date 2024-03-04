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
#include "sdf_grid.h"
#include "neuralCore/siren.h"
#include "neural_SDF/neural_sdf.h"
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

      if (positions.size() < 3*batch_size)
        positions = std::vector<float>(3*batch_size,0);
      if (distances.size() < batch_size)
        distances = std::vector<float>(batch_size,0);
      if (dparams.size() < gen_params.size()*batch_size)
        dparams = std::vector<float>(gen_params.size()*batch_size,0);
      if (dpositions.size() < 3*batch_size)
        dpositions = std::vector<float>(3*batch_size,0);

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

      if (positions.size() < 3*batch_size)
        positions = std::vector<float>(3*batch_size,0);
      if (distances.size() < batch_size)
        distances = std::vector<float>(batch_size,0);

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
    unsigned batch_size = 256;
    std::vector<float> positions, distances, dparams, dpositions;
  };

  class FieldSdfLoss : public UPGOptimizableFunction
  {
  private:
    const PointCloudReference &reference;
    unsigned batch_size = 512;
    std::vector<float> positions, distances, dparams, dpositions, ref_distances;
  public:
    FieldSdfLoss(const PointCloudReference &_reference, const Block &optimization_blk):
    reference(_reference)
    {
      
    }

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      //set parameters to ProceduralSdf
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);
      
      //check of reference and input data are ok
      assert(reference.d_points.size() > 0);  
      assert(reference.d_points.size() == reference.d_distances.size());    
      assert(gen_params.size() == out_grad.size());

      //resize temporary vectors if needed
      if (positions.size() < 3*batch_size)
        positions = std::vector<float>(3*batch_size,0);
      if (distances.size() < batch_size)
        distances = std::vector<float>(batch_size,0);
      if (ref_distances.size() < batch_size)
        ref_distances = std::vector<float>(batch_size,0);
      if (dparams.size() < gen_params.size()*batch_size)
        dparams = std::vector<float>(gen_params.size()*batch_size,0);
      if (dpositions.size() < 3*batch_size)
        dpositions = std::vector<float>(3*batch_size,0);

      unsigned param_count = gen_params.size();
      std::vector<double> out_grad_d(param_count, 0);
      double loss = 0;
      double reg_loss = 0;

      //main step - minimize SDF values on surface
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.d_points.size();
        positions[3*i+0] = reference.d_points[index].x;
        positions[3*i+1] = reference.d_points[index].y;
        positions[3*i+2] = reference.d_points[index].z;
        ref_distances[i] = reference.d_distances[index];
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), dparams.data(), dpositions.data());

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = distances[i] - ref_distances[i];
        for (int j=0;j<param_count;j++)
          out_grad_d[j] += 2*d*dparams[i*param_count + j];
        loss += d*d;
      }

      for (int i=0;i<gen_params.size();i++)
        out_grad[i] = out_grad_d[i]/batch_size;

      /*
      debug("params1 [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",params[i]);
      debug("]\n");
      debug("grad1 [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",out_grad[i]);
      debug("]\n");  
      */
     
      return loss/batch_size;
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      //set parameters to ProceduralSdf
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);
      
      //check of reference and input data are ok
      assert(reference.d_points.size() > 0);  
      assert(reference.d_points.size() == reference.d_distances.size());    

      //resize temporary vectors if needed
      if (positions.size() < 3*batch_size)
        positions = std::vector<float>(3*batch_size,0);
      if (distances.size() < batch_size)
        distances = std::vector<float>(batch_size,0);
      if (ref_distances.size() < batch_size)
        ref_distances = std::vector<float>(batch_size,0);

      unsigned param_count = gen_params.size();
      std::vector<double> out_grad_d(param_count, 0);
      double loss = 0;
      double reg_loss = 0;

      //main step - minimize SDF values on surface
      for (unsigned i=0;i<batch_size;i++)
      {
        unsigned index = rand() % reference.d_points.size();
        positions[3*i+0] = reference.d_points[index].x;
        positions[3*i+1] = reference.d_points[index].y;
        positions[3*i+2] = reference.d_points[index].z;
        ref_distances[i] = reference.d_distances[index];
      }

      sdf.get_distance_batch(batch_size, positions.data(), distances.data(), nullptr, nullptr);

      for (unsigned i=0;i<batch_size;i++)
      {
        double d = distances[i] - ref_distances[i];
        loss += d*d;
      }

      return loss/batch_size;
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
    logerr("FieldSdfLoss should not be used for constructive reconstruction!");
    return 0;
  }
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

  std::vector<UPGReconstructionResult> simple_field_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                        const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    FieldSdfLoss opt_func(reference, *step_blk);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer(step_blk->get_string("optimizer_name", "adam"), &opt_func, step_blk,
                                                            prev_step_res[0], prev_step_res[0].structure);
    optimizer->optimize();
    return optimizer->get_best_results();
  }

  std::vector<UPGReconstructionResult> sdf_grid_reconstruction(PointCloudReference &reference, uint16_t grid_node_type)
  {
    UPGParametersRaw params;
      ProceduralSdf g_sdf({{grid_node_type}});
      unsigned voxels_cnt = g_sdf.all_params.size();
      unsigned N = round(pow(voxels_cnt, 1.0/3));
      unsigned p_offset = 3+1;//move and scale params
      params.p.resize(voxels_cnt+p_offset, 0.0f);
      std::vector<float> counts(voxels_cnt, 0);
      float *voxels = params.p.data() + p_offset;

      //bbox must be square, as grid is expected the have same
      //number of voxels by each dimension and we can only
      //scale SDF by one number in every dimention.
      //we calculate scale and shift for grid transform
      AABB bbox = get_point_cloud_bbox(reference.d_points);
      float max_s = 0.5f*MAX(bbox.size().x, MAX(bbox.size().y, bbox.size().z));
      glm::vec3 bbox_center = 0.5f*(bbox.min_pos + bbox.max_pos);
      bbox = AABB(bbox_center - glm::vec3(max_s), bbox_center + glm::vec3(max_s));
      params.p[0] = bbox_center.x; //shift
      params.p[1] = bbox_center.y;
      params.p[2] = bbox_center.z;
      params.p[3] = max_s;         //scale
      
      for (int i=0;i<reference.d_points.size();i++)
      {
        glm::vec3 v = (float)N*(reference.d_points[i] - bbox.min_pos)/bbox.size();
        glm::ivec3 vi = v;
        glm::vec3 vf = v - glm::vec3(vi);
        if (v.x > 0 && v.x<N && v.y > 0 && v.y<N && v.z > 0 && v.z<N)
        {
          unsigned id = vi.x*N*N + vi.y*N + vi.z;
          if (vi.x > 0 && vi.x<N-1 && vi.y > 0 && vi.y<N-1 && vi.z > 0 && vi.z<N-1)
          {
            voxels[id]         += (1-vf.x)*(1-vf.y)*(1-vf.z)*reference.d_distances[i];
            counts[id]         += (1-vf.x)*(1-vf.y)*(1-vf.z);

            voxels[id+1]       += (1-vf.x)*(1-vf.y)*vf.z*reference.d_distances[i];
            counts[id+1]       += (1-vf.x)*(1-vf.y)*vf.z;

            voxels[id+N]       += (1-vf.x)*vf.y*(1-vf.z)*reference.d_distances[i];
            counts[id+N]       += (1-vf.x)*vf.y*(1-vf.z);
            
            voxels[id+N+1]     += (1-vf.x)*vf.y*vf.z*reference.d_distances[i];
            counts[id+N+1]     += (1-vf.x)*vf.y*vf.z;

            voxels[id+N*N]     += vf.x*(1-vf.y)*(1-vf.z)*reference.d_distances[i];
            counts[id+N*N]     += vf.x*(1-vf.y)*(1-vf.z);
            
            voxels[id+N*N+1]   += vf.x*(1-vf.y)*vf.z*reference.d_distances[i];
            counts[id+N*N+1]   += vf.x*(1-vf.y)*vf.z;

            voxels[id+N*N+N]   += vf.x*vf.y*(1-vf.z)*reference.d_distances[i];
            counts[id+N*N+N]   += vf.x*vf.y*(1-vf.z);
            
            voxels[id+N*N+N+1] += vf.x*vf.y*vf.z*reference.d_distances[i];
            counts[id+N*N+N+1] += vf.x*vf.y*vf.z;
          }
          else
          {
            voxels[id] += reference.d_distances[i];
            counts[id] += 1;
          }
        }
      }

      for (int i=0;i<voxels_cnt;i++)
      {
        voxels[i] = counts[i] > 0 ? voxels[i]/counts[i] : 1.0f/N;
      }
    /*
      GridSdfNode grid(N, bbox);
      grid.set_param_span(voxels, 0);
      std::vector<float> ddist_dpos(3, 0.0f);
      std::vector<float> ddist_dp(params_cnt, 0.0f);
      std::vector<float> x_grad(params_cnt, 0.0f);
      std::vector<float> V(params_cnt, 0.0f);
      std::vector<float> S(params_cnt, 0.0f);
      std::vector<float> &X = voxels;
      float alpha = 0.01;
      float beta_1 = 0.9;
      float beta_2 = 0.999;
      float eps = 1e-8;
      float loss = 0;
      bool verbose = true;

      unsigned samples = 10000;
      unsigned iterations = 50;
      
      for (int iter=0; iter< iterations; iter++)
      {
        std::fill_n(x_grad.begin(), x_grad.size(), 0.0f);
        double total_loss = 0.0;
        for (int i=0;i<samples;i++)
        {
          unsigned index = rand() % reference.d_points.size();
          glm::vec3 p = reference.d_points[index];
          //printf("%f %f %f -- %f\n",p.x,p.y,p.z,reference.d_distances[index]);
          //p = {urand(),urand(),urand()};
          //p = bbox.size()*p + bbox.min_pos;
          std::fill_n(ddist_dp.begin(), ddist_dp.size(), 0.0);
          float d = grid.get_distance(p, ddist_dp, ddist_dpos) - reference.d_distances[i];
          total_loss += d*d;
          for (int j=0;j<ddist_dp.size();j++)
            x_grad[j] += 2*d*ddist_dp[j] / samples;
        }
        loss = total_loss/samples;

        const float reg_q = 1e-6f;
        float reg_loss = 0;
        //regularization
        for (int i=1;i<N-1;i++)
        {
          for (int j=1;j<N-1;j++)
          {
            for (int k=1;k<N-1;k++)
            {
              unsigned bid = i*N*N + j*N + k;
              reg_loss += reg_q*(abs(X[bid+N*N+N+1] - X[bid-N*N-N-1]) + 
                                 abs(X[bid+N*N+N-1] - X[bid-N*N-N+1]) +
                                 abs(X[bid+N*N-N+1] - X[bid-N*N+N-1]) + 
                                 abs(X[bid+N*N-N-1] - X[bid-N*N+N+1]));
              x_grad[bid+N*N+N+1] += reg_q*glm::sign(X[bid+N*N+N+1] - X[bid-N*N-N-1]);
              x_grad[bid-N*N-N-1] -= reg_q*glm::sign(X[bid+N*N+N+1] - X[bid-N*N-N-1]);
              x_grad[bid+N*N+N-1] += reg_q*glm::sign(X[bid+N*N+N-1] - X[bid-N*N-N+1]);
              x_grad[bid-N*N-N+1] -= reg_q*glm::sign(X[bid+N*N+N-1] - X[bid-N*N-N+1]);
              x_grad[bid+N*N-N+1] += reg_q*glm::sign(X[bid+N*N-N+1] - X[bid-N*N+N-1]);
              x_grad[bid-N*N+N-1] -= reg_q*glm::sign(X[bid+N*N-N+1] - X[bid-N*N+N-1]);
              x_grad[bid+N*N-N-1] += reg_q*glm::sign(X[bid+N*N-N-1] - X[bid-N*N+N+1]);
              x_grad[bid-N*N+N+1] -= reg_q*glm::sign(X[bid+N*N-N-1] - X[bid-N*N+N+1]);
            }            
          }          
        }

        printf("reg loss %f\n",reg_loss);

        //debug("grad  [");
        for (int i=0;i<params_cnt;i++)
        {
          //debug("%f ",x_grad[i]);
          float g = x_grad[i];
          V[i] = beta_1 * V[i] + (1-beta_1)*g;
          float Vh = V[i] / (1 - pow(beta_1, iter+1)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
          float Sh = S[i] / (1 - pow(beta_2, iter+1)); 
          X[i] -= alpha*Vh/(sqrt(Sh) + eps);
        }
        //debug("]\n");
        if ((iter % 5 == 0) && verbose)
          debug("Adam iter %3d  val = %.8f\n", iter, loss);
      }
    */

    UPGReconstructionResult res;
    res.structure = {{SdfNodeType::MOVE, SdfNodeType::SCALE, grid_node_type}};
    res.parameters = params;
    res.loss_optimizer = 0.0f;

    //debug("params ");
    //for (auto &p : res.parameters.p)
    //  debug("%f ", p);
    //debug("\n");
/*
    int steps = 15;
    g_sdf.set_parameters(res.parameters.p);
    for (int i=0;i<steps;i++)
    {
      CameraSettings camera;
      camera.origin = glm::vec3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(g_sdf, camera, 512, 512, 16, SDFRenderMode::LAMBERT);
      engine::textureManager->save_png(t, "dataset_image_grid_"+std::to_string(i));
    }
*/
    return {res};
  }

  std::vector<UPGReconstructionResult> neural_sdf_reconstruction(PointCloudReference &reference)
  {

      //nn::TensorProcessor::init("GPU");
      nn::Siren network(nn::Siren::Type::SDF, 2, 32);
      unsigned points = reference.d_distances.size();


      //Calculate scale and shift for neural_sdf transform
      //and move all points to the unit cube
      AABB bbox = get_point_cloud_bbox(reference.d_points);
      float max_s = 0.5f*MAX(bbox.size().x, MAX(bbox.size().y, bbox.size().z));
      glm::vec3 bbox_center = 0.5f*(bbox.min_pos + bbox.max_pos);
      bbox = AABB(bbox_center - glm::vec3(max_s), bbox_center + glm::vec3(max_s));

      std::vector<float> positions(3*reference.d_distances.size(), 0);
      for (int i=0;i<points;i++)
      {
        positions[3*i+0] = (reference.d_points[i].x - bbox_center.x)/max_s;
        positions[3*i+1] = (reference.d_points[i].y - bbox_center.y)/max_s;
        positions[3*i+2] = (reference.d_points[i].z - bbox_center.z)/max_s;
      }

      network.train(positions, reference.d_distances, 512, 10000, false);

      UPGReconstructionResult res;
      res.structure = {{SdfNodeType::MOVE, SdfNodeType::SCALE, SdfNodeType::NEURAL_SMALL}};
      res.parameters.p = std::vector<float>(network.get_weights().size() + 3 + 1);
      res.parameters.p[0] = bbox_center.x;//shift
      res.parameters.p[1] = bbox_center.y;
      res.parameters.p[2] = bbox_center.z;
      res.parameters.p[3] = max_s;        //scale   
      res.parameters.p.insert(res.parameters.p.begin()+4, network.get_weights().begin(), network.get_weights().end());   
      res.loss_optimizer = 0.0f;

      /*int steps = 15;
      ProceduralSdf g_sdf(res.structure);
      g_sdf.set_parameters(res.parameters.p);
      for (int i=0;i<steps;i++)
      {
        CameraSettings camera;
        camera.origin = glm::vec3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
        camera.target = glm::vec3(0,0,0);
        camera.up = glm::vec3(0,1,0);
        Texture t = render_sdf(g_sdf, camera, 256,256,1);
        engine::textureManager->save_png(t, "image_neural_"+std::to_string(i));
      }*/

      return {res};
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
    //TODO: transform points to fit into unit cube
    
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
      else if (step_blk->get_bool("field"))
        opt_res = simple_field_reconstruction_step(step_blk, reference, opt_res);
      else if (step_blk->get_bool("grid"))
        opt_res = sdf_grid_reconstruction(reference, step_blk->get_int("grid_node_type", (int)SdfNodeType::GRID_32));
      else if (step_blk->get_bool("neural"))
        opt_res = neural_sdf_reconstruction(reference);      
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