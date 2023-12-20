#pragma once

#include "dmesh.h"
#include "camera.h"
#include <memory>

namespace diff_render
{
// /**
// \brief API to ray-scene intersection on CPU
// */
// struct CRT_Hit 
// {
//   float    t;         ///< intersection distance from ray origin to object
//   uint32_t primId; 
//   uint32_t instId;
//   uint32_t geomId;    ///< use 4 most significant bits for geometry type; thay are zero for triangles 
//   float    coords[4]; ///< custom intersection data; for triangles coords[0] and coords[1] stores baricentric coords (u,v)
// };

struct SurfaceInfo
{
  uint32_t primId; 
  uint32_t instId;
  uint32_t geomId; 
  float    t;       ///<! dist from origin ray to surface
  float    u;       ///<! first triangle baricentric 
  float    v;       ///<! second triangle baricentric 
};

struct IRayTracer
{
  IRayTracer(){}
  virtual ~IRayTracer(){}

  virtual void        Init(const Scene* pScene) = 0;
  virtual void        SetCamera(const CamInfo& cam)   = 0;
  virtual float3      GetCameraPos() const            = 0;
  virtual SurfaceInfo CastSingleRay(float x, float y, float3* outPos = nullptr, float3* outDir = nullptr) = 0;
  virtual SurfaceInfo GetNearestHit(float3 rayPos, float3 rayDir, float tNear = 0.0f, float tFar = 1e9f) = 0;
};

::std::shared_ptr<IRayTracer> MakeRayTracer3D(const char* className);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float3 EyeRayDirNormalized(float x, float y, LiteMath::float4x4 a_mViewProjInv)
{
  float4 pos = float4(2.0f*x - 1.0f, -2.0f*y + 1.0f, 0.0f, 1.0f );
  pos = a_mViewProjInv * pos;
  pos /= pos.w;
  return normalize(to_float3(pos));
}

static inline float3 mul3x3(LiteMath::float4x4 m, float3 v)
{ 
  return to_float3(m*to_float4(v, 0.0f));
}

static inline float3 mul4x3(LiteMath::float4x4 m, float3 v)
{
  return to_float3(m*to_float4(v, 1.0f));
}

static inline void transform_ray3f(LiteMath::float4x4 a_mWorldViewInv, float3* ray_pos, float3* ray_dir) 
{
  float3 pos  = mul4x3(a_mWorldViewInv, (*ray_pos));
  float3 pos2 = mul4x3(a_mWorldViewInv, ((*ray_pos) + 100.0f*(*ray_dir)));

  float3 diff = pos2 - pos;

  (*ray_pos)  = pos;
  (*ray_dir)  = normalize(diff);
}
}