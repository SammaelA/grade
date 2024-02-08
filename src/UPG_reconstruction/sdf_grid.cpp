#include "sdf_grid.h"

namespace upg
{
  float GridSdfNode::get_distance(const glm::vec3 &pos, std::vector<float> *ddist_dp,
                                  std::vector<float> *ddist_dpos) const
  {
    glm::vec3 vox_f = grid_size_f*glm::min(glm::max((pos-bbox.min_pos)/bbox_size, 0.0f), 1.0f-1e-5f);
    glm::uvec3 vox_u = vox_f;
    glm::vec3 dp = vox_f - glm::vec3(vox_u);

    float res = 0.0;
    
    #define id(i,j,k) (vox_u.z + i)*grid_size*grid_size + (vox_u.y + j)*grid_size + (vox_u.x + k)
    if (vox_u.x<grid_size-1 && vox_u.y<grid_size-1 && vox_u.z<grid_size-1)
    {
      res += p[id(0,0,0)]*(1 - dp.x)*(1 - dp.y)*(1 - dp.z);
      res += p[id(0,0,1)]*(1 - dp.x)*(1 - dp.y)*(dp.z);
      res += p[id(0,1,0)]*(1 - dp.x)*(dp.y)*(1 - dp.z);
      res += p[id(0,1,1)]*(1 - dp.x)*(dp.y)*(dp.z);
      res += p[id(1,0,0)]*(dp.x)*(1 - dp.y)*(1 - dp.z);
      res += p[id(1,0,1)]*(dp.x)*(1 - dp.y)*(dp.z);
      res += p[id(1,1,0)]*(dp.x)*(dp.y)*(1 - dp.z);
      res += p[id(1,1,1)]*(dp.x)*(dp.y)*(dp.z);
    }
    else
      res += p[id(0,0,0)];


    return res;
  }
  std::vector<ParametersDescription::Param> GridSdfNode::get_parameters_block(AABB scene_bbox) const
  {
    return {};
  }
  void GridSdfNode::set_voxel(const glm::uvec3 &vox, float distance)
  {
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
  void GridSdfNode::set_voxel(const glm::vec3 &pos, float distance)
  {
    glm::uvec3 vox = grid_size_f*glm::min(glm::max((pos-bbox.min_pos)/bbox_size, 0.0f), 1.0f);
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
  float GridSdfNode::sample_bilinear(const glm::vec3 &pos, std::vector<float> *ddist_dp,
                                     std::vector<float> *ddist_dpos) const
  {
    glm::vec3 vox_f = grid_size_f*glm::min(glm::max((pos-bbox.min_pos)/bbox_size, 0.0f), 1.0f);
    glm::uvec3 vox_u = vox_f;
    glm::vec3 dp = vox_f - glm::vec3(vox_u);

    float res = 0.0;
    
    #define id(i,j,k) (vox_u.z + i)*grid_size*grid_size + (vox_u.y + j)*grid_size + (vox_u.x + k)
    res += p[id(0,0,0)]*(1 - dp.x)*(1 - dp.y)*(1 - dp.z);
    res += p[id(0,0,1)]*(1 - dp.x)*(1 - dp.y)*(dp.z);
    res += p[id(0,1,0)]*(1 - dp.x)*(dp.y)*(1 - dp.z);
    res += p[id(0,1,1)]*(1 - dp.x)*(dp.y)*(dp.z);
    res += p[id(1,0,0)]*(dp.x)*(1 - dp.y)*(1 - dp.z);
    res += p[id(1,0,1)]*(dp.x)*(1 - dp.y)*(dp.z);
    res += p[id(1,1,0)]*(dp.x)*(dp.y)*(1 - dp.z);
    res += p[id(1,1,1)]*(dp.x)*(dp.y)*(dp.z);

    return res;
  }
}