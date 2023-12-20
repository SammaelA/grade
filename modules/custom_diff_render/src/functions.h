#pragma once
#include "raytrace.h"
namespace diff_render
{
static inline void VertexShader(const CamInfo& u, float vx, float vy, float vz, 
                                float output[2])
{
  const float W    =   vx * u.mWVP[3] + vy * u.mWVP[7] + vz * u.mWVP[11] + u.mWVP[15]; 
  const float xNDC =  (vx * u.mWVP[0] + vy * u.mWVP[4] + vz * u.mWVP[ 8] + u.mWVP[12])/W;
  const float yNDC = -(vx * u.mWVP[1] + vy * u.mWVP[5] + vz * u.mWVP[ 9] + u.mWVP[13])/W;
  output[0] = (xNDC*0.5f + 0.5f)*u.width;
  output[1] = (yNDC*0.5f + 0.5f)*u.height; 
}

void VS_X_grad(float V[3], const CamInfo &data, float _d_V[3]);
void VS_Y_grad(float V[3], const CamInfo &data, float _d_V[3]);

void BarU_grad(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3], 
               float* _d_A, float* _d_B, float* _d_C);

void BarV_grad(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3], 
               float* _d_A, float* _d_B, float* _d_C);

}