#include "sdf_grid.h"

namespace upg
{
  void GridSdfNode::get_distance_batch(unsigned     batch_size,
                                       float *const positions,    
                                       float *      distances,
                                       float *      ddist_dparams,
                                       float *      ddist_dpos,
                               std::vector<float> & stack,
                                       unsigned     stack_head) const
  {
    if (ddist_dparams)
    {
      for (int i=0;i<batch_size;i++)
        distances[i] = get_distance(glm::vec3(positions[3*i+0], positions[3*i+1], positions[3*i+2]), 
                                    std::span<float>(ddist_dparams + batch_size*p_offset + i*param_cnt(), param_cnt()),
                                    std::span<float>(ddist_dpos + 3*i, 3));
    }
    else
    {
      for (int i=0;i<batch_size;i++)
        distances[i] = get_distance(glm::vec3(positions[3*i+0], positions[3*i+1], positions[3*i+2]), std::span<float>(), std::span<float>()); 
    }
  }

  float GridSdfNode::get_distance(const glm::vec3 &pos, std::span<float> ddist_dp, 
                                  std::span<float> ddist_dpos) const
  {
    //bbox for grid is a unit cube
    glm::vec3 vox_f = grid_size_f*glm::min(glm::max((pos-glm::vec3(-1,-1,-1))/bbox_size, 0.0f), 1.0f-1e-5f);
    glm::uvec3 vox_u = vox_f;
    glm::vec3 dp = vox_f - glm::vec3(vox_u);

    float res = 0.0;
    
    #define id(i,j,k) (vox_u.z + k)*grid_size*grid_size + (vox_u.y + j)*grid_size + (vox_u.x + i)
    #define C(i,j,k) {\
      float qx = (1 - dp.x + i*(2*dp.x-1));\
      float qy = (1 - dp.y + j*(2*dp.y-1));\
      float qz = (1 - dp.z + k*(2*dp.z-1));\
      res += p[id(i,j,k)]*qx*qy*qz;\
      if (ddist_dp.data()){\
        ddist_dp[id(i,j,k)] += qx*qy*qz;\
        ddist_dpos[0] += p[id(i,j,k)]*qy*qz * (2*i-1);\
        ddist_dpos[1] += p[id(i,j,k)]*qx*qz * (2*j-1);\
        ddist_dpos[2] += p[id(i,j,k)]*qx*qy * (2*k-1);\
      }\
    }
    if (vox_u.x<grid_size-1 && vox_u.y<grid_size-1 && vox_u.z<grid_size-1)
    {
      C(0,0,0);
      C(0,0,1);
      C(0,1,0);
      C(0,1,1);
      C(1,0,0);
      C(1,0,1);
      C(1,1,0);
      C(1,1,1);
    }
    else
    {  
      res += p[id(0,0,0)];
      if (ddist_dp.data())
        ddist_dp[id(0,0,0)] += 1;

      //no dependency on position on edges
      //(*ddist_dpos)[0] += 0;
      //(*ddist_dpos)[1] += 0;
      //(*ddist_dpos)[2] += 0;
    }

    //printf("dist %f\n",res);

    return res;
  }
  std::vector<ParametersDescription::Param> GridSdfNode::get_parameters_block(AABB scene_bbox) const
  {
    std::vector<ParametersDescription::Param> params;
    params.push_back({1,-1.0f,1.0f, ParameterType::ARRAY, "data", grid_size*grid_size*grid_size});
    return params;
  }
  void GridSdfNode::set_voxel(const glm::uvec3 &vox, float distance)
  {
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
  void GridSdfNode::set_voxel(const glm::vec3 &pos, float distance)
  {
    //bbox for grid is a unit cube
    glm::uvec3 vox = grid_size_f*glm::min(glm::max((pos-glm::vec3(-1,-1,-1))/bbox_size, 0.0f), 1.0f);
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
}