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

    
    //  Calculate cube border to take points near pos to calculate interpolation result 
    //  X_left, X_right, Y_left, Y_Right, Z_left, Z_right
    
    std::vector<uint32_t> sample_bbox {0, 0, 0, 0, 0, 0};

    int sample_cube_size = 4;
    
    //  Axis X
    if ((float)vox_u.x - (float)sample_cube_size / 2 < 0)
    {
      sample_bbox[0] = 0;
      sample_bbox[1] = sample_cube_size - 1;
    }
    else if ((float)vox_u.x + (float)sample_cube_size / 2 > 31)
    {
      sample_bbox[1] = 31;
      sample_bbox[0] = 31 - sample_cube_size + 1;
    }
    else
    {
      sample_bbox[0] = (float)vox_u.x - (float)sample_cube_size / 2 + 1;
      sample_bbox[1] = (float)vox_u.x + (float)sample_cube_size / 2;
    }

    //  Axis Y
    if ((float)vox_u.y - (float)sample_cube_size / 2 < 0)
    {
      sample_bbox[2] = 0;
      sample_bbox[3] = sample_cube_size - 1;
    }
    else if ((float)vox_u.y + (float)sample_cube_size / 2 > 31)
    {
      sample_bbox[3] = 31;
      sample_bbox[2] = 31 - sample_cube_size + 1;
    }
    else
    {
      sample_bbox[2] = (float)vox_u.y - (float)sample_cube_size / 2 + 1;
      sample_bbox[3] = (float)vox_u.y + (float)sample_cube_size / 2;
    }

    //  Axis Z
    if ((float)vox_u.z - (float)sample_cube_size / 2 < 0)
    {
      sample_bbox[4] = 0;
      sample_bbox[5] = sample_cube_size - 1;
    }
    else if ((float)vox_u.z + (float)sample_cube_size / 2 > 31)
    {
      sample_bbox[5] = 31;
      sample_bbox[4] = 31 - sample_cube_size + 1;
    }
    else
    {
      sample_bbox[4] = (float)vox_u.z - (float)sample_cube_size / 2 + 1;
      sample_bbox[5] = (float)vox_u.z + (float)sample_cube_size / 2;
    }

    //  Save points to find interpolation coefs
    std::vector<LiteMath::float3> X;
    std::vector<float> b;

    for (int x = sample_bbox[0]; x <= sample_bbox[1]; x++)
    {
      for (int y = sample_bbox[2]; y <= sample_bbox[3]; y++)
      {
        for (int z = sample_bbox[4]; z <= sample_bbox[5]; z++)
        {
          int index = (z * grid_size + y) * grid_size + x;

          X.push_back(LiteMath::float3(x, y, z) * bbox_size / grid_size_f + LiteMath::float3(-1, -1, -1));
          b.push_back(p[index]);
        }
      }
    }

    auto A = interpolation::create_A(X);
    std::vector<float> Q(64 * 64, 0), R(64 * 64, 0);
    interpolation::householder_qr(A, 64, Q, R);

    std::vector<float> G = interpolation::mul_qr(Q, R, 64);
    std::vector<float> coefs = interpolation::calc_qr_coefs(Q, R, b);

    res = interpolation::perform_interpolation(coefs, pos);

    // std::cout << res << " " << p[id(0, 0, 0)] << std::endl;
    // std::cout << interpolation::matrix_norm(A, interpolation::mul_qr(Q, R, 64)) << std::endl;

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