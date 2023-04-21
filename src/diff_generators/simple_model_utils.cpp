#include <vector>
#include <string>
#include "diff_geometry_generation.h"
#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/transform.hpp>
#include "simple_model_utils.h"

namespace dgen
{
  void shift(std::vector<float> &model, glm::vec3 sh)
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      model[i*FLOAT_PER_VERTEX] += sh.x;
      model[i*FLOAT_PER_VERTEX+1] += sh.y;
      model[i*FLOAT_PER_VERTEX+2] += sh.z;
    }
  }

  void scale(std::vector<float> &model, glm::vec3 scale)
  {
    glm::mat4 scale_mat = glm::scale(glm::mat4(1.0f),scale); 
    transform(model, scale_mat);
  }

  void transform(std::vector<float> &model, glm::mat4 transform_mat)
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      glm::vec4 p = glm::vec4(model[i*FLOAT_PER_VERTEX], model[i*FLOAT_PER_VERTEX+1], model[i*FLOAT_PER_VERTEX+2], 1);
      glm::vec4 n4 = glm::vec4(model[i*FLOAT_PER_VERTEX+3], model[i*FLOAT_PER_VERTEX+4], model[i*FLOAT_PER_VERTEX+5], 0);

      glm::mat4 n_mat = glm::transpose(glm::inverse(transform_mat));

      p = transform_mat*p;
      n4 = n_mat*n4;
      glm::vec3 n = glm::normalize(glm::vec3(n4.x, n4.y, n4.z));

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
    glm::vec3 min_pos(1e18,1e18,1e18);
    glm::vec3 max_pos(-1e18,-1e18,-1e18);

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

  glm::vec3 get_center_of_mass_vertices(const std::vector<float> &model)
  {
    glm::vec3 CoM(0,0,0);
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      CoM += glm::vec3(model[i*FLOAT_PER_VERTEX], model[i*FLOAT_PER_VERTEX+1], model[i*FLOAT_PER_VERTEX+2]);
    }
    return CoM/MAX(1.0f, model.size() / FLOAT_PER_VERTEX);
  }

  void normalize_model(std::vector<float> &model)
  {
    AABB bbox = get_bbox(model);
    glm::vec3 sizes = bbox.max_pos - bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    shift(model, -0.5f*(bbox.max_pos + bbox.min_pos));
    scale(model, glm::vec3(1/max_size));
  }
}