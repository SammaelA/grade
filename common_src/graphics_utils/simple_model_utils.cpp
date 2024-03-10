#include <vector>
#include <string>
//#include "diff_geometry_generation.h"
#include "common_utils/matrix_transform.h"
#include "simple_model_utils.h"
#include "common_utils/blk.h"

namespace visualizer
{
  static constexpr int FLOAT_PER_VERTEX = 8;

  void shift(std::vector<float> &model, float3 sh)
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      model[i*FLOAT_PER_VERTEX] += sh.x;
      model[i*FLOAT_PER_VERTEX+1] += sh.y;
      model[i*FLOAT_PER_VERTEX+2] += sh.z;
    }
  }

  void scale(std::vector<float> &model, float3 scale)
  {
    float4x4 scale_mat = LiteMath::scale(float4x4(),scale); 
    transform(model, scale_mat);
  }

  void transform(std::vector<float> &model, float4x4 transform_mat)
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      float4 p = float4(model[i*FLOAT_PER_VERTEX], model[i*FLOAT_PER_VERTEX+1], model[i*FLOAT_PER_VERTEX+2], 1);
      float4 n4 = float4(model[i*FLOAT_PER_VERTEX+3], model[i*FLOAT_PER_VERTEX+4], model[i*FLOAT_PER_VERTEX+5], 0);

      float4x4 n_mat = LiteMath::transpose(LiteMath::inverse4x4(transform_mat));

      p = transform_mat*p;
      n4 = n_mat*n4;
      float3 n = normalize(float3(n4.x, n4.y, n4.z));

      model[i*FLOAT_PER_VERTEX] = p.x;
      model[i*FLOAT_PER_VERTEX+1] = p.y;
      model[i*FLOAT_PER_VERTEX+2] = p.z;

      model[i*FLOAT_PER_VERTEX+3] = n.x;
      model[i*FLOAT_PER_VERTEX+4] = n.y;
      model[i*FLOAT_PER_VERTEX+5] = n.z;
    }
  }
  AABB get_bbox(const std::vector<float> &model)
  {
    float3 min_pos(1e18,1e18,1e18);
    float3 max_pos(-1e18,-1e18,-1e18);

    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      min_pos.x = MIN(min_pos.x, model[i*FLOAT_PER_VERTEX]);
      min_pos.y = MIN(min_pos.y, model[i*FLOAT_PER_VERTEX+1]);
      min_pos.z = MIN(min_pos.z, model[i*FLOAT_PER_VERTEX+2]);
      
      max_pos.x = MAX(max_pos.x, model[i*FLOAT_PER_VERTEX]);
      max_pos.y = MAX(max_pos.y, model[i*FLOAT_PER_VERTEX+1]);
      max_pos.z = MAX(max_pos.z, model[i*FLOAT_PER_VERTEX+2]);
    }

    return AABB(min_pos, max_pos);
  }

  float3 get_center_of_mass_vertices(const std::vector<float> &model)
  {
    float3 CoM(0,0,0);
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      CoM += float3(model[i*FLOAT_PER_VERTEX], model[i*FLOAT_PER_VERTEX+1], model[i*FLOAT_PER_VERTEX+2]);
    }
    return CoM/MAX(1.0f, model.size() / FLOAT_PER_VERTEX);
  }

  void normalize_model(std::vector<float> &model)
  {
    AABB bbox = get_bbox(model);
    float3 sizes = bbox.max_pos - bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    shift(model, -0.5f*(bbox.max_pos + bbox.min_pos));
    scale(model, float3(1/max_size));
  }

  void save_camera_settings(const CameraSettings &camera, Block &blk)
  {
    blk.set_vec3("camera.origin", camera.origin);
    blk.set_vec3("camera.target", camera.target);
    blk.set_vec3("camera.up", camera.up);
    blk.set_double("camera.z_near", camera.z_near);
    blk.set_double("camera.z_far", camera.z_far);
    blk.set_double("camera.fov_rad", camera.fov_rad);
    blk.set_bool("camera.orthographic", camera.orthographic);
  }

  CameraSettings load_camera_settings(const Block &blk)
  {
    CameraSettings camera;
    
    camera.origin = blk.get_vec3("camera.origin", camera.origin);
    camera.target = blk.get_vec3("camera.target", camera.target);
    camera.up = blk.get_vec3("camera.up", camera.up);
    camera.z_near = blk.get_double("camera.z_near", camera.z_near);
    camera.z_far = blk.get_double("camera.z_far", camera.z_far);
    camera.fov_rad = blk.get_double("camera.fov_rad", camera.fov_rad);
    camera.orthographic = blk.get_bool("camera.orthographic", camera.orthographic);

    return camera;
  }
}