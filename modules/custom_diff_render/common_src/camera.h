#pragma once
#include "LiteMath.h"

namespace diff_render
{
struct CamInfo
{
  CamInfo() = default;
  CamInfo(float3 _origin, float3 _target, float3 _up, float w, float h, float _fov_rad = M_PI_4, float _z_near = 0.01, float _z_far = 100):
  origin(_origin),
  target(_target),
  up(_up),
  width(w),
  height(h),
  fov_rad(_fov_rad),
  zNear(_z_near),
  zFar(_z_far)
  {
    mProj = LiteMath::perspectiveMatrix(fov_rad/LiteMath::DEG_TO_RAD, width/height, zNear, zFar);
    mWorldView = LiteMath::lookAt(origin, target, up);
    commit();
  }
  float3 origin;
  float3 target;
  float3 up;
  float fov_rad = M_PI_4;
  float zNear = 0.01;
  float zFar = 100;
  LiteMath::float4x4 mWorldView;
  LiteMath::float4x4 mProj;

  float mWVP[16]; // WorlViewProject := (mProj*(mView*mWorld))
  float width;
  float height;
  
  /**
  \brief make all needed internal computations, prepare for rendering
  */
  void commit()
  {
    LiteMath::float4x4 mTransform = mProj*mWorldView;
    memcpy(mWVP, (float*)&mTransform, 16*sizeof(float));
  }
};
} // namespace diff_render
