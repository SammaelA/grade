#include "upg.h"
#include "sdf_node.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "common_utils/bbox.h"
#include "sdf_rendering.h"
#include "graphics_utils/render_point_cloud.h"
#include "common_utils/distribution.h"

namespace upg
{

  struct PointCloudReference
  {
    std::vector<glm::vec3> points;
    std::vector<glm::vec3> outside_points;

    std::vector<glm::vec3> d_points;
    std::vector<float> d_distances;

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
      sdf_to_point_cloud_with_dist(sdf, 2*points, &(reference.d_points), &(reference.d_distances));

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

  float normal_pdf(float x, float x0, float sigma)
  {
    return (1 / sqrt(2 * PI * sigma * sigma)) * exp(-SQR(x - x0) / (2 * sigma * sigma));
  }

  // estimate sum(f(0), ..., f(max_index)) by sampling only part of these values
  // function expect all f(i) >= 0 and increases number of samples if average f(i) is small
  // returns partial sum and number of samples
  std::pair<int, double> sum_with_adaptive_batching(std::function<double(int)> get_value, 
                                                             int max_index, int base_batch_size,
                                                             int max_batch_size,
                                                             double sensitivity)
  {
    int samples = 0;
    double partial_sum = 0.0;
    if (base_batch_size < max_index)
    {
      while (samples < max_batch_size && partial_sum <= sensitivity)
      {
        for (int b = 0; b < base_batch_size; b++)
        {
          int index = rand() % max_index;
          partial_sum += get_value(index);
        }
        samples += base_batch_size;
        sensitivity /= 10;
        if (sensitivity == 0)
          break;
      }
    }
    else
    {
      samples = max_index;
      for (int i = 0; i < max_index; i++)
        partial_sum += get_value(i);
    }

    return {samples, partial_sum};
  }

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
        d = glm::sign(d)*std::min(0.03, abs(d));
        for (int i=0;i<gen_params.size();i++)
          out_grad_d[i] += 2*d*cur_grad[i];
        return d*d;
      }, reference.points.size(), 256, 1024, 10);

      //regularization - penalty for outside points with sdf < 0
      auto p2 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.outside_points[index];
        cur_grad.clear();
        double d = sdf.get_distance(p, &cur_grad, &dpos_dparams);
        if (d < 0)
        {
          for (int i=0;i<gen_params.size();i++)
            out_grad_d[i] += (d > 0 ? 2 : 4)*d*cur_grad[i];
          return d*d;
        }
        else
          return 0;
      }, reference.outside_points.size(), 128, 1024, 0.1);

      /*
      debug("params [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",params[i]);
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
        d = glm::sign(d)*std::min(0.03, abs(d));
        return d*d;
      }, reference.points.size(), 256, 1024, 10);

      //regularization - penalty for outside points with sdf < 0
      auto p2 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &p = reference.outside_points[index];
        double d = sdf.get_distance(p);
        if (d < 0)
          return d*d;
        else
          return 0;
      }, reference.outside_points.size(), 128, 1024, 0.1);

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

    virtual float estimate_positioning_quality(const UPGStructure &structure,
                                               const UPGPart &part, std::span<const float> parameters,
                                               float border_sigma,
                                               float inner_point_penalty) const override
  {
    UPGStructure part_structure;
    part_structure.s = std::vector<uint16_t>(structure.s.begin() + part.s_range.first, structure.s.begin() + part.s_range.second);
    std::span<const float> part_parameters(parameters.data() + part.p_range.first, parameters.data() + part.p_range.second);

    SdfGenInstance gen(part_structure);
    ProceduralSdf sdf = gen.generate(part_parameters);

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

    logerr("part [%d %d][%d %d]", (int)part.s_range.first, (int)part.s_range.second, part.p_range.first, part.p_range.second);
    logerr("quality %f (%f %f)", (float)quality, (float)q1, (float)q2);

    return quality;
  }
  private:
    const PointCloudReference &reference;
  };

  class ConstructiveSdfCompare : public UPGOptimizableFunction
  {
  private:
    std::vector<glm::vec3> surface_points;
    std::vector<glm::vec3> outer_points;
    float beta;
    float p;
    float beta_p;
    float beta_p_inv;
    float inner_fine;

  public:
    ConstructiveSdfCompare(const std::vector<glm::vec3> &_surface_points,
                           const std::vector<glm::vec3> &_outer_points,
                           float _beta = 0.1, float _p = 2, float _inner_fine = 1)
    {
      surface_points = _surface_points;
      outer_points = _outer_points;
      beta = _beta;
      p = _p;
      beta_p = pow(beta, p);
      beta_p_inv = 1/beta_p;
      inner_fine = _inner_fine;
    }

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(surface_points.size() > 0);
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
        const glm::vec3 &point = surface_points[index];
        cur_grad.clear();
        double dist = sdf.get_distance(point, &cur_grad, &dpos_dparams);

        //loss function, value in range (0, 1]
        double dist_inv = 1/MAX(abs(dist), 1e-12);
        double dist_p = pow(abs(dist),p);
        double dist_p_inv = 1/dist_p;
        double loss = abs(dist) >= beta ? (1 - 0.5*beta_p*dist_p_inv) : (0.5*beta_p_inv*dist_p);
        double dloss_ddist = abs(dist) >= beta ? (0.5*beta_p*p*(dist_p_inv*dist_inv)*glm::sign(dist)) : (0.5*beta_p_inv*(dist_p*dist_inv)*glm::sign(dist));

        //if (dist < -beta)
        //{
        //  loss += inner_fine*(-dist-beta);
        //  dloss_ddist += -inner_fine;
        //}

        for (int i=0;i<gen_params.size();i++)
          out_grad_d[i] += dloss_ddist*cur_grad[i];
        return loss;
      }, surface_points.size(), 512, 512, 1);

      //regularization - penalty for outside points with sdf < 0
      std::pair<int, double> p2 = {0,0};
      if (outer_points.size() > 0)
      {
        p2 = sum_with_adaptive_batching([&](int index) -> double 
        {
          const glm::vec3 &p = outer_points[index];
          cur_grad.clear();
          double d = sdf.get_distance(p, &cur_grad, &dpos_dparams);
          if (d < 0)
          {
            for (int i=0;i<gen_params.size();i++)
              out_grad_d[i] += 2*inner_fine*d*cur_grad[i];
            return inner_fine*d*d;
          }
          else
            return 0;
        }, outer_points.size(), 512, 512, 1);
      }
      /*
      debug("params [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",params[i]);
      debug("]\n");
      debug("grad [");
      for (int i=0;i<gen_params.size();i++)
        debug("%f ",out_grad[i]);
      debug("]\n");
      */
     //logerr("P1 P2 %f %f\n", p1.second/p1.first, p2.second);
      
      for (int i=0;i<gen_params.size();i++)
        out_grad[i] = out_grad_d[i]/(p1.first+p2.first);
      
      return p1.second/p1.first + p2.second/MAX(1, p2.first);
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(surface_points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);

      //main step - minimize SDF values on surface
      auto p1 = sum_with_adaptive_batching([&](int index) -> double 
      {
        const glm::vec3 &point = surface_points[index];
        double dist = sdf.get_distance(point);

        //loss function, value in range (0, 1]
        double dist_inv = 1/MAX(abs(dist), 1e-12);
        double dist_p = pow(abs(dist),p);
        double dist_p_inv = 1/dist_p;
        double loss = abs(dist) >= beta ? (1 - 0.5*beta_p*dist_p_inv) : (0.5*beta_p_inv*dist_p);
        //if (dist < -beta)
        //  loss += inner_fine*(-dist-beta);
        return loss;
      }, surface_points.size(), 512, 512, 1);

      //regularization - penalty for outside points with sdf < 0
      std::pair<int, double> p2 = {0,0};
      if (outer_points.size() > 0)
      {
        p2 = sum_with_adaptive_batching([&](int index) -> double 
        {
          const glm::vec3 &p = outer_points[index];
          double d = sdf.get_distance(p);
          if (d < 0)
            return inner_fine*d*d;
          else
            return 0;
        }, outer_points.size(), 512, 512, 1);
      }

      return p1.second/p1.first + p2.second/MAX(1, p2.first);
    }
    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) override
    {
      return ((SdfGenInstance*)gen)->desc;
    }
    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const override
    {
      return std::make_shared<SdfGenInstance>(structure);
    }
  };

struct ReferencePointsGrid
{
public:
  ReferencePointsGrid() = default;
  ReferencePointsGrid(const std::vector<glm::vec3> &_points, const std::vector<float> &_distances,
                      unsigned cells_count);

  struct Cell
  {
    unsigned offset;//in arrays below
    unsigned count;
    unsigned active_count;
    AABB bbox;
  };

  inline unsigned v_to_i(glm::uvec3 v)
  {
    return v.x*grid_size.y*grid_size.z + v.y*grid_size.z + v.z;
  }

  inline glm::uvec3 get_cell_id(glm::vec3 point)
  {
    glm::vec3 idf = (point-grid_box.min_pos)/cell_size;
    return glm::uvec3(CLAMP(floor(idf.x), 0, grid_size.x-1),
                      CLAMP(floor(idf.y), 0, grid_size.y-1),
                      CLAMP(floor(idf.z), 0, grid_size.z-1));
  }

  //same size for all three arrays
  std::vector<glm::vec3> points;
  std::vector<float> distances;
  std::vector<bool> is_active;

  unsigned cells_count;
  glm::uvec3 grid_size;
  glm::vec3 cell_size;
  AABB grid_box;
  std::vector<Cell> cells;
};

ReferencePointsGrid::ReferencePointsGrid(const std::vector<glm::vec3> &_points, const std::vector<float> &_distances,
                                         unsigned expected_cells_count)
{
  assert(_points.size() > 0);
  assert(_distances.size() == _points.size());
  assert(expected_cells_count < 10'000'000); //too large grid may be result of some bug

  grid_box = get_point_cloud_bbox(_points);
  glm::vec3 box_size = grid_box.max_pos - grid_box.min_pos;
  float min_size = MIN(box_size.x, MIN(box_size.y, box_size.z));
  float k_flt = pow(expected_cells_count/((box_size.x/min_size)*(box_size.y/min_size)*(box_size.z/min_size)), 1/3.0);
  if (min_size < 1e-4 || k_flt < 0.5)
  {
    logerr("ReferencePointsGrid: give point cloud has size (%f %f %f), extremely small size in one dimension can lead to issues!", 
           box_size.x, box_size.y, box_size.z);
  }

  grid_size = glm::uvec3(ceil(k_flt*(box_size.x/min_size)), ceil(k_flt*(box_size.y/min_size)), ceil(k_flt*(box_size.z/min_size)));
  cell_size = box_size/glm::vec3(grid_size);
  cells_count = grid_size.x*grid_size.y*grid_size.z;
  logerr("ReferencePointsGrid: %u cells (%u %u %u)", cells_count, grid_size.x, grid_size.y, grid_size.z);
  cells.resize(cells_count);

  std::vector<std::vector<unsigned>> point_lists;
  point_lists.resize(cells_count);
  for (unsigned i=0;i<_points.size();i++)
  {
    unsigned cid = v_to_i(get_cell_id(_points[i]));
    point_lists[cid].push_back(i);
  }

  points = std::vector<glm::vec3>(_points.size());
  distances = std::vector<float>(_points.size());
  is_active = std::vector<bool>(_points.size(), true);
  //points = _points;
  //distances = _distances;
  //return;

  unsigned offset = 0;
  for (unsigned i=0;i<cells_count;i++)
  {
    for (int j=0;j<point_lists[i].size();j++)
    {
      points[offset+j] = _points[point_lists[i][j]];
      distances[offset+j] = _distances[point_lists[i][j]];
    }
    cells[i].offset = offset;
    cells[i].count = point_lists[i].size();
    cells[i].active_count = point_lists[i].size();

    glm::uvec3 cell_id(i/(grid_size.y*grid_size.z), (i/grid_size.z)%grid_size.y, i%grid_size.z);
    glm::vec3 pos = grid_box.min_pos + glm::vec3(cell_id)*cell_size;
    cells[i].bbox = AABB(pos, pos+cell_size);

    //logerr("cell %u (off %u). %u points. bbox [%f %f %f]-[%f %f %f]", i, offset, cells[i].count, 
    //       cells[i].bbox.min_pos.x, cells[i].bbox.min_pos.y, cells[i].bbox.min_pos.z,
    //       cells[i].bbox.max_pos.x, cells[i].bbox.max_pos.y, cells[i].bbox.max_pos.z);
    offset += point_lists[i].size();
  }
}

class FieldSdfCompare : public UPGOptimizableFunction
  {
  public:
    ReferencePointsGrid grid;
    std::vector<glm::vec3> a_points;
    std::vector<float> a_distances;

    std::vector<float> ddist_dparams, ddist_dpos;
    std::vector<double> ddist_dparams_sum;

    float beta;
    float p;
    float beta_p;
    float beta_p_inv;
    float reg_q;
    int samples;

  public:
    bool verbose = false;

    FieldSdfCompare(const std::vector<glm::vec3> &_points,
                    const std::vector<float> &_distances,
                    unsigned grid_cells_count = 16*16*16,
                    unsigned max_parameters = 1024):
//                    float _beta = 0.1, float _p = 2, float _reg_q = 1e6)
    grid(_points, _distances, grid_cells_count)
    {
      a_points = _points;
      a_distances = _distances;
      ddist_dparams.resize(max_parameters);
      ddist_dparams_sum.resize(max_parameters);
      ddist_dpos = {0,0,0};
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
      std::vector<glm::vec3> active_points;
      std::vector<float> active_distances;
      for (int i=0;i<a_points.size();i++)
      {
        if (abs(sdf.get_distance(a_points[i]) - a_distances[i]) >= threshold)
        {
          active_points.push_back(a_points[i]);
          active_distances.push_back(a_distances[i]);
        }
      }
      logerr("removed %d/%d points", (int)(a_points.size()-active_points.size()), (int)a_points.size());
      a_points = active_points;
      a_distances = active_distances;

      unsigned active_left = 0;
      unsigned active_cells = 0;
      for (auto &c : grid.cells)
      {
        if (c.active_count == 0)
          continue;
        for (unsigned i=c.offset;i<c.offset+c.count;i++)
        {
          if (grid.is_active[i] && abs(sdf.get_distance(grid.points[i]) - grid.distances[i]) < threshold)
          {
            grid.is_active[i] = false;
            c.active_count--;
          }
        }
        active_left += c.active_count;
        active_cells += (c.active_count>0);
      }
      logerr("active cells left %u points %u", active_cells, active_left);
    }

    template <bool is_differentiable>
    double main_loss_func(const ProceduralSdf &sdf, const glm::vec3 &point, float target_dist)
    {
      ddist_dparams.clear();
      double dist = sdf.get_distance(point, &ddist_dparams, &ddist_dpos) - target_dist;

      // loss function, value in range (0, 1]
      double dist_p = pow(abs(dist), p);
      double dist_p_inv = 1 / dist_p;
      double loss = abs(dist) >= beta ? (-0.5 * beta_p * dist_p_inv) : (0.5 * beta_p_inv * dist_p - 1);

      if constexpr (is_differentiable)
      {
        double dist_inv = 1 / MAX(abs(dist), 1e-12);
        double dloss_ddist = abs(dist) >= beta ? (0.5 * beta_p * p * (dist_p_inv * dist_inv) * glm::sign(dist)) : (0.5 * beta_p_inv * (dist_p * dist_inv) * glm::sign(dist));
        for (int i = 0; i < ddist_dparams_sum.size(); i++)
          ddist_dparams_sum[i] += dloss_ddist * ddist_dparams[i];
      }
      return loss;
    }

    template <bool is_differentiable>
    double reg_loss_func(const ProceduralSdf &sdf, const glm::vec3 &point, float target_dist)
    {
      ddist_dparams.clear();
      double dist = sdf.get_distance(point, &ddist_dparams, &ddist_dpos) - target_dist;
      dist = glm::sign(dist)*MAX(0.05, abs(dist));
      
      if (dist < 0)
      {
        double loss = dist*dist;
        if constexpr (is_differentiable)
        {
          for (int i = 0; i < ddist_dparams_sum.size(); i++)
            ddist_dparams_sum[i] += 2*reg_q*dist*ddist_dparams[i];
        }
        return loss;
      }
      else
        return 0;
    }

    template <bool is_differentiable>
    float sample_random(const ProceduralSdf &sdf, int samples, std::span<float> out_grad = {})
    {
      if constexpr (is_differentiable)
        ddist_dparams_sum = std::vector<double>(out_grad.size(), 0);
      
      double total_loss = 0;
      int Na=1;
      if (samples < 0.75*a_points.size())
      {
        Na = samples;
        for (int i=0;i<samples;i++)
        {
          int index = rand() % a_points.size();
          total_loss += main_loss_func<is_differentiable>(sdf, a_points[index], a_distances[index]);
        }
      }
      else
      {
        Na = a_points.size();
        for (int i=0;i<a_points.size();i++)
          total_loss += main_loss_func<is_differentiable>(sdf, a_points[i], a_distances[i]);
      }

      if constexpr (is_differentiable)
      {
        for (int i=0;i<out_grad.size();i++)
          out_grad[i] = ddist_dparams_sum[i]/Na;
        ddist_dparams_sum = std::vector<double>(out_grad.size(), 0);
      }

      double total_loss_r = 0;
      int Ns = 1;
      if (samples < 0.75*grid.points.size())
      {
        Ns = samples;
        for (int i=0;i<samples;i++)
        {
          int index = rand() % grid.points.size();
          total_loss_r += reg_loss_func<is_differentiable>(sdf, grid.points[index], grid.distances[index]);
        }
      }
      else
      {
        Ns = grid.points.size();
        for (int i=0;i<grid.points.size();i++)
          total_loss_r += reg_loss_func<is_differentiable>(sdf, grid.points[i], grid.distances[i]);
      }
      if constexpr (is_differentiable)
      {
        for (int i=0;i<out_grad.size();i++)
          out_grad[i] += ddist_dparams_sum[i]/Ns;
      }
      
      if (verbose)
        logerr("loss %lg %lg", total_loss, total_loss_r);

      return total_loss/Na + total_loss_r/Ns;
    }

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

    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd,
                           const OptParams &params, std::span<float> out_grad) override
    {
      assert(a_points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);
      assert(gen_params.size() == out_grad.size());

      return sample_random<true>(sdf, samples, out_grad);
    }

    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) override
    {
      assert(a_points.size() > 0);
      std::vector<float> gen_params = opt_params_to_gen_params(params, pd);
      ProceduralSdf sdf = ((SdfGenInstance*)gen)->generate(gen_params);
      return sample_random<false>(sdf, samples);
    }

    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) override
    {
      return ((SdfGenInstance*)gen)->desc;
    }

    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const override
    {
      return std::make_shared<SdfGenInstance>(structure);
    }
  };

  std::shared_ptr<UPGOptimizer> get_optimizer(const std::string &optimizer_name,
                                              UPGOptimizableFunction *opt_func,
                                              Block *step_blk,
                                              UPGReconstructionResult start_params = UPGReconstructionResult(),
                                              UPGStructure structure = UPGStructure())
  {
    std::shared_ptr<UPGOptimizer> optimizer;
    if (optimizer_name == "adam")
      optimizer = get_optimizer_adam(opt_func, *step_blk, start_params);
    else if (optimizer_name == "memetic")
      optimizer = get_optimizer_memetic(opt_func, *step_blk, structure);
    else if (optimizer_name == "CHC")
      optimizer = get_optimizer_CHC(opt_func, *step_blk, structure);
    else if (optimizer_name == "particle_swarm")
      optimizer = get_optimizer_particle_swarm(opt_func, *step_blk, structure);
    else if (optimizer_name == "CC")
      optimizer = get_optimizer_CC(opt_func, *step_blk, structure);
    else if (optimizer_name == "DE")
      optimizer = get_optimizer_differentiable_evolution(opt_func, *step_blk, structure);
    else if (optimizer_name == "iterative_fitting")
      optimizer = get_optimizer_iterative_fitting(opt_func, *step_blk, structure);
    return optimizer;
  }

  std::vector<UPGReconstructionResult> simple_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                  const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    SdfRenderAndCompare opt_func(reference, *step_blk);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer(step_blk->get_string("optimizer_name", "adam"), &opt_func, step_blk,
                                                            prev_step_res[0], prev_step_res[0].structure);
    optimizer->optimize();
    return optimizer->get_best_results();
  }

  UPGReconstructionResult merge_results(const std::vector<UPGReconstructionResult> &part_results)
  {
    std::vector<uint16_t> structure_inv;
    std::vector<float> params_inv;
    uint32_t num = 1;
    for (auto &pr : part_results)
    {
      for (int i=pr.structure.s.size()-1;i>=0;i--)
        structure_inv.push_back(pr.structure.s[i]);
      for (int i=pr.parameters.p.size()-1;i>=0;i--)
        params_inv.push_back(pr.parameters.p[i]);      
      int32_t S = 1;
      while (num>= S && ((num & S) == 0))
      {
        structure_inv.push_back(3); //merge node id
        S = S << 1;
      }
      num++;
    }
    int32_t S = 1;
    while (num>= S && ((num & S) == 0))
    {
      structure_inv.push_back(3); //merge node id
      S = S << 1;
    }

    UPGReconstructionResult res;
    for (int i=0;i<structure_inv.size();i++)
      res.structure.s.push_back(structure_inv[structure_inv.size()-i-1]);
    for (int i=0;i<params_inv.size();i++)
      res.parameters.p.push_back(params_inv[params_inv.size()-i-1]);
    debug("res structure {");
    for (auto &s : res.structure.s)
      debug("%d ",(int)s);
    debug("}\n");
    return res;
  }

  void reassign_point_types(std::vector<glm::vec3> &surface_points/*inout*/, std::vector<glm::vec3> &outer_points/*inout*/, 
                            ProceduralSdf &sdf, float threshold = 0.001f)
  {
    std::vector<glm::vec3> new_surface_points;
    std::vector<glm::vec3> new_outer_points = outer_points;

    float  beta = 0.01f;
    float  P = 2.0f;
    float  beta_p = pow(beta, P);
    float  beta_p_inv = 1/beta_p;

    int reassigned = 0;
    int target = 0;
    float in_beta = 0.0f;
    float out_beta = 0.0f;
    for (auto &p : surface_points)
    {
      if (abs(0.5-length(p-glm::vec3(0.6,0.6,0.6))) < 1e-5) 
        target++;
      float dist = sdf.get_distance(p);
        double dist_inv = 1/MAX(abs(dist), 1e-12);
        double dist_p = pow(abs(dist),P);
        double dist_p_inv = 1/dist_p;
      float loss = abs(dist) >= beta ? (1 - 0.5*beta_p*dist_p_inv) : (0.5*beta_p_inv*dist_p);
      if (abs(dist) < beta)
        in_beta += 1-loss;
      else
        out_beta+= 1-loss;
      if (abs(sdf.get_distance(p)) < threshold)
      {
        new_outer_points.push_back(p);
        reassigned++; 
      }
      else
        new_surface_points.push_back(p);
    }

    surface_points = new_surface_points;
    outer_points = new_outer_points;
    logerr("inbeta out beta %f %f", in_beta, out_beta);
    logerr("reassigned %d/%d points %d/%d surface points left", reassigned, target, (int)surface_points.size(), (int)(surface_points.size()+outer_points.size()));
  }

  void remove_inactive_points(std::vector<glm::vec3> &points, std::vector<float> &distances, ProceduralSdf &sdf, float threshold = 0.001f)
  {
    std::vector<glm::vec3> active_points;
    std::vector<float> active_distances;
    for (int i=0;i<points.size();i++)
    {
      if (abs(sdf.get_distance(points[i]) - distances[i]) >= threshold)
      {
        active_points.push_back(points[i]);
        active_distances.push_back(distances[i]);
      }
    }
    logerr("removed %d/%d points", (int)(points.size()-active_points.size()), (int)points.size());
    points = active_points;
    distances = active_distances;
  }

  std::vector<UPGReconstructionResult> field_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                 const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    std::vector<UPGReconstructionResult> part_results;
    std::vector<glm::vec3> active_points = reference.d_points;
    std::vector<float> active_distances = reference.d_distances;

    //step_blk->set_int("iterations", 1000);
    //step_blk->set_bool("verbose", true);
    //step_blk->set_double("learning_rate", 0.005f);
    //step_blk->set_string("optimizer_name", "adam");
    Block fine_opt_blk;
    fine_opt_blk.copy(step_blk);
    fine_opt_blk.set_string("optimizer_name", "adam");
    fine_opt_blk.set_int("iterations", 1000);
    fine_opt_blk.set_bool("verbose", true);
    fine_opt_blk.set_double("learning_rate", 0.0001f);

    auto parts = get_sdf_parts(prev_step_res[0].structure);
    FieldSdfCompare opt_func(reference.d_points, reference.d_distances, 8*8*8);

    for (int i=0; i<parts.size(); i++)
    {      
      srand(time(NULL));
      UPGStructure part_structure;
      for (int j=parts[i].s_range.first; j<parts[i].s_range.second; j++)
      {
        part_structure.s.push_back(prev_step_res[0].structure.s[j]);
      }

      opt_func.set_hyperparameters(0.01f, 1.5f, 2048, 50000);
      opt_func.verbose = false;
      std::shared_ptr<UPGOptimizer> optimizer = get_optimizer(step_blk->get_string("optimizer_name", "adam"), &opt_func, step_blk,
                                                              UPGReconstructionResult(), part_structure);
      optimizer->optimize();
      auto partial_result = optimizer->get_best_results()[0];
      logerr("res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                  partial_result.parameters.p[2], partial_result.parameters.p[3]);

      opt_func.set_hyperparameters(0.002f, 1.5f, reference.d_points.size(), 50000);
      opt_func.verbose = true;
      std::shared_ptr<UPGOptimizer> fine_optimizer = get_optimizer("adam", &opt_func, &fine_opt_blk,
                                                                   partial_result, partial_result.structure);
      fine_optimizer->optimize();
      partial_result = fine_optimizer->get_best_results()[0];
      logerr("tuned res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                        partial_result.parameters.p[2], partial_result.parameters.p[3]);

      part_results.push_back(partial_result);

      SdfGenInstance gen(partial_result.structure);
      ProceduralSdf sdf = gen.generate(partial_result.parameters.p);
      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 4);
      engine::textureManager->save_png(t, "partial_sdf_"+std::to_string(i));
      opt_func.remove_inactive_points(sdf, 0.001f);
    }
    return {merge_results(part_results)};

  }

  UPGReconstructionResult set_initial_parameters(const UPGStructure &structure, FieldSdfCompare &func, glm::vec3 initial_pos)
  {
    auto gen = func.get_generator(structure);
    auto pd = func.get_full_parameters_description(gen.get());

    glm::ivec3 pos_ids(-1,-1,-1);
    std::vector<glm::vec2> borders;
    for (const auto &p : pd.get_block_params())
    {
      for (const auto &param_info : p.second.p)
      {
        if (param_info.type != ParameterType::CONST)
          borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
        if (param_info.name=="move_x")
        {
          assert(pos_ids.x == -1);
          pos_ids.x = borders.size()-1;
        }
        else if (param_info.name=="move_y")
        {
          assert(pos_ids.y == -1);
          pos_ids.y = borders.size()-1;
        }
        else if (param_info.name=="move_z")
        {
          assert(pos_ids.z == -1);
          pos_ids.z = borders.size()-1;
        }
      }
    }
    assert(pos_ids.x>=0 && pos_ids.y>=0 && pos_ids.z>=0);

    UPGReconstructionResult res;
    res.structure = structure;
    res.parameters.p.resize(borders.size());
    float range = 0.5;

    for (int i=0;i<borders.size();i++)
    {
      float sz = borders[i].y - borders[i].x;
      if (borders[i].x < 0)
        res.parameters.p[i] = urand(borders[i].x + 0.5*(1-range)*sz, borders[i].y - 0.5*(1-range)*sz);
      else
        res.parameters.p[i] = urand(borders[i].x, borders[i].x + 0.25*sz);
    }

    res.parameters.p[pos_ids.x] = initial_pos.x;
    res.parameters.p[pos_ids.y] = initial_pos.y;
    res.parameters.p[pos_ids.z] = initial_pos.z;

    return res;
  }

  std::vector<UPGReconstructionResult> field_reconstruction_step2(Block *step_blk, PointCloudReference &reference,
                                                                  const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    std::vector<UPGReconstructionResult> part_results;
    int tries = 15;

    Block opt_blk;
    opt_blk.copy(step_blk);
    opt_blk.set_string("optimizer_name", "adam");
    opt_blk.set_int("iterations", 1000);
    opt_blk.set_bool("verbose", false);
    opt_blk.set_double("learning_rate", 0.005f);

    Block fine_opt_blk;
    fine_opt_blk.copy(step_blk);
    fine_opt_blk.set_string("optimizer_name", "adam");
    fine_opt_blk.set_int("iterations", 200);
    fine_opt_blk.set_bool("verbose", true);
    fine_opt_blk.set_double("learning_rate", 0.0005f);

    auto parts = get_sdf_parts(prev_step_res[0].structure);
    FieldSdfCompare opt_func(reference.d_points, reference.d_distances, 8*8*8);

    for (int i=0; i<parts.size(); i++)
    {      
      srand(time(NULL));
      UPGStructure part_structure;
      for (int j=parts[i].s_range.first; j<parts[i].s_range.second; j++)
        part_structure.s.push_back(prev_step_res[0].structure.s[j]);
      opt_func.set_hyperparameters(0.01f, 2.0f, 2048, 100);
      opt_func.verbose = false;

      UPGReconstructionResult best_part_result;
      best_part_result.loss_optimizer = 1e9;

      for (int try_n=0;try_n<tries; try_n++)
      {
        AABB bbox = opt_func.grid.grid_box;
        glm::vec3 best_pos;
        float best_dist = 1e6;
        for (int k=0;k<15000;k++)
        {
          unsigned index = rand() % opt_func.grid.points.size();
          if (opt_func.grid.is_active[index] && opt_func.grid.distances[index] < best_dist)
          {
            best_dist = opt_func.grid.distances[index];
            best_pos = opt_func.grid.points[index];
          }
        }

        logerr("start pos %f %f %f", best_pos.x, best_pos.y, best_pos.z);
        auto init_params = set_initial_parameters(part_structure, opt_func, best_pos);
        logerr("start val = %f %f %f %f",init_params.parameters.p[0], init_params.parameters.p[1],
                                         init_params.parameters.p[2], init_params.parameters.p[3]);
        std::shared_ptr<UPGOptimizer> optimizer = get_optimizer("adam", &opt_func, &opt_blk, init_params, part_structure);
        optimizer->optimize();
        auto partial_result = optimizer->get_best_results()[0];
        logerr("res %f = %f %f %f %f",partial_result.loss_optimizer, partial_result.parameters.p[0], partial_result.parameters.p[1],
                                    partial_result.parameters.p[2], partial_result.parameters.p[3]);
        if (partial_result.loss_optimizer < best_part_result.loss_optimizer)
          best_part_result = partial_result;
      }
      logerr("best res = %f %f %f %f", best_part_result.parameters.p[0], best_part_result.parameters.p[1],
                                       best_part_result.parameters.p[2], best_part_result.parameters.p[3]);

      opt_func.set_hyperparameters(0.002f, 2.0f, reference.d_points.size(), 100);
      opt_func.verbose = true;
      std::shared_ptr<UPGOptimizer> fine_optimizer = get_optimizer("adam", &opt_func, &fine_opt_blk,
                                                                   best_part_result, best_part_result.structure);
      fine_optimizer->optimize();
      best_part_result = fine_optimizer->get_best_results()[0];
      logerr("tuned res = %f %f %f %f", best_part_result.parameters.p[0], best_part_result.parameters.p[1],
                                        best_part_result.parameters.p[2], best_part_result.parameters.p[3]);

      part_results.push_back(best_part_result);

      SdfGenInstance gen(best_part_result.structure);
      ProceduralSdf sdf = gen.generate(best_part_result.parameters.p);
      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 4);
      engine::textureManager->save_png(t, "partial_sdf_"+std::to_string(i));
      opt_func.remove_inactive_points(sdf, 0.002f);
    }
    return {merge_results(part_results)};

  }

  std::vector<UPGReconstructionResult> constructive_reconstruction_step(Block *step_blk, PointCloudReference &reference,
                                                                        const std::vector<UPGReconstructionResult> &prev_step_res)
  {
    std::vector<UPGReconstructionResult> part_results;
    std::vector<glm::vec3> surface_points = reference.points;
    std::vector<glm::vec3> outer_points = reference.outside_points;

    //step_blk->set_int("iterations", 1000);
    //step_blk->set_bool("verbose", true);
    //step_blk->set_double("learning_rate", 0.005f);
    //step_blk->set_string("optimizer_name", "adam");
    Block fine_opt_blk;
    fine_opt_blk.copy(step_blk);
    fine_opt_blk.set_string("optimizer_name", "adam");
    fine_opt_blk.set_int("iterations", 1000);
    fine_opt_blk.set_bool("verbose", false);
    fine_opt_blk.set_double("learning_rate", 0.001f);

    auto parts = get_sdf_parts(prev_step_res[0].structure);

    for (int i=0; i<parts.size(); i++)
    {      
      srand(time(NULL));
      UPGStructure part_structure;
      for (int j=parts[i].s_range.first; j<parts[i].s_range.second; j++)
      {
        logerr("struct %d",  prev_step_res[0].structure.s[j]);
        part_structure.s.push_back(prev_step_res[0].structure.s[j]);
      }
      //part_structure.s = {3,2,1,2,1};
      UPGReconstructionResult pr;
      //pr.structure = part_structure;
      //pr.parameters.p = {0.5,0.7,0.56,0.3};

      ConstructiveSdfCompare opt_func(surface_points, outer_points, 0.01f, 1.5f, 1e12);
      std::shared_ptr<UPGOptimizer> optimizer = get_optimizer(step_blk->get_string("optimizer_name", "adam"), &opt_func, step_blk,
                                                              pr, part_structure);
      optimizer->optimize();
      auto partial_result = optimizer->get_best_results()[0];
      logerr("res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                  partial_result.parameters.p[2], partial_result.parameters.p[3]);

      ConstructiveSdfCompare fine_opt_func(surface_points, outer_points, 0.002f, 1.5f, 1e12);
      std::shared_ptr<UPGOptimizer> fine_optimizer = get_optimizer("adam", &fine_opt_func, &fine_opt_blk,
                                                                   partial_result, partial_result.structure);
      fine_optimizer->optimize();
      partial_result = fine_optimizer->get_best_results()[0];
      logerr("tuned res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                        partial_result.parameters.p[2], partial_result.parameters.p[3]);

      part_results.push_back(partial_result);

      SdfGenInstance gen(partial_result.structure);
      ProceduralSdf sdf = gen.generate(partial_result.parameters.p);
      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 4);
      engine::textureManager->save_png(t, "partial_sdf_"+std::to_string(i));
      reassign_point_types(surface_points, outer_points, sdf);
    }
    return {merge_results(part_results)};

  }

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
      if (step_blk->get_bool("constructive_reconstruction"))
        opt_res = field_reconstruction_step2(step_blk, reference, opt_res);
      else
        opt_res = simple_reconstruction_step(step_blk, reference, opt_res);
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
        result.quality_ir = get_sdf_similarity_MSE(reference_sdf, sdf);
      }
      if (reference.is_synthetic && res_blk->get_bool("check_image_quality"))
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