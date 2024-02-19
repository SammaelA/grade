#pragma once
#include "upg.h"
#include "sdf_node.h"
#include "optimization.h"

namespace upg
{
  class FieldSdfCompare : public UPGOptimizableFunction
  {
  public:
    const std::vector<glm::vec3> &points;
    const std::vector<float> &distances;
    //ReferencePointsGrid grid;
    std::vector<glm::vec3> a_points;
    std::vector<float> a_distances;
    int active_points_count = 0;

    std::vector<float> ddist_dparams, ddist_dpos;
    std::vector<double> ddist_dparams_sum;
    std::vector<float> distances_batch, a_distances_batch;
    std::vector<glm::vec3> batch;

    float beta;
    float p;
    float beta_p;
    float beta_p_inv;
    float reg_q;
    int samples;

    unsigned long seed = 0ul;

  public:
    bool verbose = false;

    FieldSdfCompare(const std::vector<glm::vec3> &_points,
                    const std::vector<float> &_distances,
                    unsigned grid_cells_count = 16*16*16,
                    unsigned max_parameters = 1024):
//                    float _beta = 0.1, float _p = 2, float _reg_q = 1e6)
    points(_points),
    distances(_distances)
//    grid(_points, _distances, grid_cells_count)
    {
      a_points = _points;
      a_distances = _distances;
      active_points_count = a_points.size();
      seed = rand();
      distances_batch.resize(_distances.size());
      a_distances_batch.resize(_distances.size());
      batch.resize(_distances.size());
    }

    //fast and multiprocessing-safe random for stochastic gradient descent
    unsigned next_index()
    {
      seed = seed * 1103515245 + 12345;
      return((unsigned)(seed/65536) % active_points_count);
    }

    void set_hyperparameters(float _beta = 0.1, float _p = 2, int _samples = 512, float _reg_q = 1e6)
    {
      beta = _beta;
      p = _p;
      beta_p = pow(beta, p);
      beta_p_inv = 1/beta_p;
      reg_q = _reg_q;
      samples = _samples;
    }

    void remove_inactive_points(ProceduralSdf &sdf, float threshold = 0.001f)
    {
      sdf.get_distance_batch(a_points.size(), (float*)a_points.data(), distances_batch.data(), nullptr, nullptr);
      std::vector<glm::vec3> active_points;
      active_points.reserve(a_points.size());
      std::vector<float> active_distances;
      active_distances.reserve(a_points.size());
      for (int i=0;i<a_points.size();i++)
      {
        if (abs(distances_batch[i] - a_distances[i]) >= threshold)
        {
          active_points.push_back(a_points[i]);
          active_distances.push_back(a_distances[i]);
        }
      }
      logerr("removed %d/%d points", (int)(a_points.size()-active_points.size()), (int)a_points.size());
      if (active_points.empty())
      {
        active_points = {a_points[0]};
        active_distances = {a_distances[0]};
      }
      a_points = active_points;
      a_distances = active_distances;
      active_points_count = a_points.size();

/*
      unsigned active_left = 0;
      unsigned active_cells = 0;
      for (auto &c : grid.cells)
      {
        if (c.active_count == 0)
          continue;
        for (unsigned i=c.offset;i<c.offset+c.count;i++)
        {
          if (grid.is_active[i] && abs(sdf.get_distance(points[i]) - distances[i]) < threshold)
          {
            grid.is_active[i] = false;
            c.active_count--;
          }
        }
        active_left += c.active_count;
        active_cells += (c.active_count>0);
      }
      logerr("active cells left %u points %u", active_cells, active_left);
*/
    }

    template <bool is_differentiable>
    double main_loss_func(float sdf_dist, float target_dist, float *sample_ddist_dparams)
    {
      float dist = sdf_dist - target_dist;
      // loss function, value in range (0, 1]
      double dist_p = pow(abs(dist), p);
      double dist_p_inv = 1 / dist_p;
      double loss = abs(dist) >= beta ? (-0.5 * beta_p * dist_p_inv) : (0.5 * beta_p_inv * dist_p - 1);

      if constexpr (is_differentiable)
      {
        double dist_inv = 1 / MAX(abs(dist), 1e-12);
        double dloss_ddist = abs(dist) >= beta ? (0.5 * beta_p * p * (dist_p_inv * dist_inv) * glm::sign(dist)) : (0.5 * beta_p_inv * (dist_p * dist_inv) * glm::sign(dist));
        for (int i = 0; i < ddist_dparams_sum.size(); i++)
          ddist_dparams_sum[i] += dloss_ddist * sample_ddist_dparams[i];
      }
      return loss;
    }

    template <bool is_differentiable>
    double reg_loss_func(float sdf_dist, float target_dist, float *sample_ddist_dparams)
    {
      float dist = sdf_dist - target_dist;
      dist = glm::sign(dist)*MAX(0.05, abs(dist));
      
      if (dist < 0)
      {
        double loss = reg_q*dist*dist;
        if constexpr (is_differentiable)
        {
          for (int i = 0; i < ddist_dparams_sum.size(); i++)
            ddist_dparams_sum[i] += 2*reg_q*dist*sample_ddist_dparams[i];
        }
        return loss;
      }
      else
        return 0;
    }

    template <bool is_differentiable>
    float sample_random(const ProceduralSdf &sdf, int samples, std::span<float> out_grad = {})
    {
      int param_cnt = out_grad.size();

      if constexpr (is_differentiable)
        ddist_dparams_sum = std::vector<double>(param_cnt, 0);
      
      double total_loss = 0;
      int Na=1;
      if (samples < 0.75*active_points_count)
      {
        Na = samples;
        for (int i=0;i<samples;i++)
        {
          int index = next_index();
          batch[i] = a_points[index];
          a_distances_batch[i] = a_distances[index];
        }
      }
      else
      {
        Na = active_points_count;
        batch = a_points;
        a_distances_batch = a_distances;
      }

      sdf.get_distance_batch(Na, (float*)batch.data(), distances_batch.data(), 
                             is_differentiable ? ddist_dparams.data() : nullptr,
                             is_differentiable ? ddist_dpos.data() : nullptr);
                             
      for (int i=0;i<Na;i++)
        total_loss += main_loss_func<is_differentiable>(distances_batch[i],a_distances_batch[i],
                                                        is_differentiable ? ddist_dparams.data() + i*param_cnt : nullptr);
      
      if constexpr (is_differentiable)
      {
        for (int i=0;i<param_cnt;i++)
          out_grad[i] = ddist_dparams_sum[i]/Na;
        ddist_dparams_sum = std::vector<double>(param_cnt, 0);
      }

      double total_loss_r = 0;
      for (int i=0;i<Na;i++)
        total_loss_r += reg_loss_func<is_differentiable>(distances_batch[i],a_distances_batch[i],
                                                         is_differentiable ? ddist_dparams.data() + i*param_cnt : nullptr);
      
      return total_loss/Na + total_loss_r/Na;
    }

/*
    template <bool is_differentiable>
    float sample_grid(const ProceduralSdf &sdf, int samples, std::span<float> out_grad = {})
    {
      if constexpr (is_differentiable)
        ddist_dparams_sum = std::vector<double>(out_grad.size(), 0);
      double total_loss = 0;
      double total_loss_r = 0;

      AABB bbox = sdf.root->get_bbox();
      glm::vec3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
      glm::vec3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
      AABB inflated_bbox = AABB(center - 1.5f*size, center + 1.5f*size);
      
      glm::uvec3 reg_st = grid.get_cell_id(inflated_bbox.min_pos);
      glm::uvec3 reg_end = min(grid.grid_size, grid.get_cell_id(inflated_bbox.max_pos) + glm::uvec3(1,1,1));
      unsigned reg_size = (reg_end.x-reg_st.x)*(reg_end.y-reg_st.y)*(reg_end.z-reg_st.z);

      //if (verbose)
      //logerr("region (%u %u %u)-(%u %u %u) %u/%u cells", reg_st.x, reg_st.y, reg_st.z, reg_end.x, reg_end.y, reg_end.z, reg_size,
      //       grid.cells_count);
      unsigned reg_points = 0;
      unsigned reg_active = 0;
      for (unsigned i=reg_st.x; i<reg_end.x; i++)
      {
        for (unsigned j=reg_st.y; j<reg_end.y; j++)
        {
          for (unsigned k=reg_st.z; k<reg_end.z; k++)
          {
            unsigned cell_id = i*grid.grid_size.y*grid.grid_size.z + j*grid.grid_size.z + k;
            reg_points += grid.cells[cell_id].count;
            reg_active += grid.cells[cell_id].active_count;
          }
        }
      }

      float sample_chance = (float)samples/reg_points;
      float active_chance = (float)samples/reg_active;

      unsigned Ns = 0;
      unsigned Na = 0;

      for (unsigned i=reg_st.x; i<reg_end.x; i++)
      {
        for (unsigned j=reg_st.y; j<reg_end.y; j++)
        {
          for (unsigned k=reg_st.z; k<reg_end.z; k++)
          {
            unsigned cell_id = i*grid.grid_size.y*grid.grid_size.z + j*grid.grid_size.z + k;
            for (unsigned index = grid.cells[cell_id].offset; index < grid.cells[cell_id].offset + grid.cells[cell_id].count; index++)
            {
              if (urand() < sample_chance)
              {
                total_loss_r += reg_loss_func<is_differentiable>(sdf, grid.points[index], grid.distances[index]);
                Ns++;
              }
              if (grid.is_active[index] && urand() < active_chance)
              {
                total_loss += main_loss_func<is_differentiable>(sdf, grid.points[index], grid.distances[index]);
                Na++;
              }
            }
          }
        }
      }

      //logerr("%lg %lg %u %u", total_loss_r, total_loss, Ns, Na);
      if constexpr (is_differentiable)
      {
        for (int i=0;i<out_grad.size();i++)
          out_grad[i] = ddist_dparams_sum[i]/(Ns + Na);
      }
      
      return total_loss/Na + total_loss_r/Ns;
      //if (verbose)
      //logerr("points %u active %u", reg_points, reg_active);

      for (int i=0;i<samples;i++)
      {
        int index = rand() % a_points.size();
        total_loss += main_loss_func<is_differentiable>(sdf, a_points[index], a_distances[index]);
      }
      for (int i=0;i<samples;i++)
      {
        int index = rand() % grid.points.size();
        total_loss_r += reg_loss_func<is_differentiable>(sdf, grid.points[index], grid.distances[index]);
      }

      if constexpr (is_differentiable)
      {
        for (int i=0;i<out_grad.size();i++)
          out_grad[i] = ddist_dparams_sum[i]/(2*samples);
      }
      
      return total_loss/samples + total_loss_r/samples;
    }
*/

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(a_points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);
      assert(gen_params.size() == out_grad.size());
      if (ddist_dparams.size() < samples*gen_params.size())
        ddist_dparams.resize(samples*gen_params.size());
      if (ddist_dparams_sum.size() < gen_params.size())
        ddist_dparams_sum.resize(gen_params.size());
      if (ddist_dpos.size() < 3*samples)
        ddist_dpos.resize(3*samples);

      return sample_random<true>(sdf, samples, out_grad);
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(a_points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf &sdf = *((ProceduralSdf*)gen);
      sdf.set_parameters(gen_params);
      return sample_random<false>(sdf, samples);
    }

    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) const override
    {
      return ((ProceduralSdf*)gen)->desc;
    }

    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const override
    {
      return std::make_shared<ProceduralSdf>(structure);
    }
  };
}