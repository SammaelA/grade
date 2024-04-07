#include "sdf_grid.h"
#include <fstream>

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
        distances[i] = get_distance(float3(positions[3*i+0], positions[3*i+1], positions[3*i+2]), 
                                    std::span<float>(ddist_dparams + batch_size*p_offset + i*param_cnt(), param_cnt()),
                                    std::span<float>(ddist_dpos + 3*i, 3));
    }
    else
    {
      for (int i=0;i<batch_size;i++)
        distances[i] = get_distance(float3(positions[3*i+0], positions[3*i+1], positions[3*i+2]), std::span<float>(), std::span<float>()); 
    }
  }

  float GridSdfNode::get_distance_svd(const float3 &pos, std::span<float> ddist_dp, 
                                  std::span<float> ddist_dpos) const
  {
    //bbox for grid is a unit cube
    float3 vox_f = grid_size_f*LiteMath::min(LiteMath::max((pos-float3(-1,-1,-1))/bbox_size, 0.0f), 1.0f-1e-5f);
    uint3 vox_u = uint3(vox_f);
    float3 dp = vox_f - float3(vox_u);

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

  float GridSdfNode::get_distance(const float3 &pos, std::span<float> ddist_dp, 
                                  std::span<float> ddist_dpos) const
  {
    //bbox for grid is a unit cube
    float3 vox_f = grid_size_f*LiteMath::min(LiteMath::max((pos-float3(-1,-1,-1))/bbox_size, 0.0f), 1.0f-1e-5f);
    uint3 vox_u = uint3(vox_f);
    float3 dp = vox_f - float3(vox_u);

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

    // if (vox_u.x<grid_size-1 && vox_u.y<grid_size-1 && vox_u.z<grid_size-1)
    // {
    //   C(0,0,0);
    //   C(0,0,1);
    //   C(0,1,0);
    //   C(0,1,1);
    //   C(1,0,0);
    //   C(1,0,1);
    //   C(1,1,0);
    //   C(1,1,1);
    // }
    // else
    // {  
    //   res += p[id(0,0,0)];
    //   if (ddist_dp.data())
    //     ddist_dp[id(0,0,0)] += 1;

    //   //no dependency on position on edges
    //   //(*ddist_dpos)[0] += 0;
    //   //(*ddist_dpos)[1] += 0;
    //   //(*ddist_dpos)[2] += 0;
    // }

    // return res;

    //printf("dist %f\n",res);


    if (vox_u.x<grid_size-1 && vox_u.y<grid_size-1 && vox_u.z<grid_size-1)
    {
      std::vector<float> b(8, 0);

      b[0] = p[id(0, 0, 0)];
      b[1] = p[id(1, 0, 0)];
      b[2] = p[id(0, 1, 0)];
      b[3] = p[id(1, 1, 0)];
      b[4] = p[id(0, 0, 1)];
      b[5] = p[id(1, 0, 1)];
      b[6] = p[id(0, 1, 1)];
      b[7] = p[id(1, 1, 1)];

      //  derrivatives
      // float dx[8], dy[8], dz[8];

      // dx[0] = b[1] - b[0];
      // dx[1] = -dx[0];
      // dx[2] = b[3] - b[2];
      // dx[3] = -dx[2];
      // dx[4] = b[5] - b[4];
      // dx[5] = -dx[4];
      // dx[6] = b[7] - b[6];
      // dx[7] = -dx[6];

      // dy[0] = b[2] - b[0];
      // dy[1] = b[3] - b[1];
      // dy[2] = -dy[0];
      // dy[3] = -dy[1];
      // dy[4] = b[6] - b[4];
      // dy[5] = b[7] - b[5];
      // dy[6] = -dy[4];
      // dy[7] = -dy[5];

      // dz[0] = b[4] - b[0];
      // dz[1] = b[5] - b[1];
      // dz[2] = b[6] - b[2];
      // dz[3] = b[7] - b[3];
      // dz[4] = -dz[0];
      // dz[5] = -dz[1];
      // dz[6] = -dz[2];
      // dz[7] = -dz[3];

      auto coefs = interpolation::calc_coefs(b);
      res = interpolation::calc_interpolation(coefs, dp);
    }
    else
    {  
      res += p[id(0,0,0)];
      // std::cout << vox_u.x << " " << vox_u.y << " " << vox_u.z << std::endl;
    }

    return res;
  }

  std::vector<ParametersDescription::Param> GridSdfNode::get_parameters_block(AABB scene_bbox) const
  {
    std::vector<ParametersDescription::Param> params;
    params.push_back({1,-1.0f,1.0f, ParameterType::ARRAY, "data", grid_size*grid_size*grid_size});
    return params;
  }
  void GridSdfNode::set_voxel(const uint3 &vox, float distance)
  {
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
  void GridSdfNode::set_voxel(const float3 &pos, float distance)
  {
    //bbox for grid is a unit cube
    uint3 vox = uint3(grid_size_f*LiteMath::min(LiteMath::max((pos-float3(-1,-1,-1))/bbox_size, 0.0f), 1.0f));
    *((float*)p.data() + vox.z*grid_size*grid_size + vox.y*grid_size + vox.x) = distance;
  }
  void GridSdfNode::save(const std::string &path) const
  {
    std::ofstream fs(path, std::ios::binary);
    unsigned count = p.size();

    fs.write((const char *)(&count), sizeof(unsigned));
    fs.write((const char *)p.data(), count * sizeof(float));
    fs.flush();
    fs.close();
  }
  void GridSdfNode::load(const std::string &path)
  {
    std::ifstream fs(path, std::ios::binary);
    unsigned count = 0;

    fs.read((char *)(&count), sizeof(unsigned));

    assert(count == p.size());

    fs.read((char *)p.data(), count * sizeof(float));
    fs.close();
  }
  void GridSdfNode::primitive_SDF_to_grid(ProceduralSdf &sdf, const AABB &bbox, float *grid, unsigned vox_size)
  {
    int offset = 0;
    for (int i = 0; i < vox_size; i++)
    {
      for (int j = 0; j < vox_size; j++)
      {
        for (int k = 0; k < vox_size; k++)
        {
          float3 p = {(k + 0.5) / vox_size, (j + 0.5) / vox_size, (i + 0.5) / vox_size};
          p = bbox.size() * p + bbox.min_pos;
          float d = sdf.get_distance(p);
          grid[offset] = d;
          offset++;
        }
      }
    }
  }
}