#pragma once
#include "common_utils/matrix_transform.h"

struct CameraSettings
{
  float3 origin, target, up;
  float fov_rad = 3.14159265f / 3;
  float z_near = 0.1f;
  float z_far = 100.0f;
  bool orthographic = false;
  
  float4x4 get_view() const
  {
    return LiteMath::lookAtRH(origin, target, up);
  }
  float4x4 get_proj(bool use_y_swap = true) const
  {
    float4x4 y_swap = to_float4x4(float4(1,0,0,0), float4(0,-1,0,0), float4(0,0,1,0),float4(0,0,0,1));
    float4x4 projection = orthographic ?
                           LiteMath::ortho(-0.5f*(z_far-z_near), 0.5f*(z_far-z_near), -0.5f*(z_far-z_near), 0.5f*(z_far-z_near), z_near, z_far) :
                           LiteMath::perspective(fov_rad, 1.0f, z_near, z_far);
    if (use_y_swap)
      return y_swap*projection;
    else
      return projection;
  }
  float4x4 get_viewProj(bool use_y_swap = true) const
  {
    return get_proj(use_y_swap)*get_view();
  }
};

struct Camera
{
    float3 pos = float3(0, 10, 10);
    float3 front = float3(0, 0, -10);
    float3 up = float3(0, 1, 0);
    float4x4 camera_mat = float4x4();
    float yaw = 0;
    float pitch = 0;
    float roll = 0;
    float4x4 &camera() { camera_mat = LiteMath::lookAtRH(pos, pos + front, up); return camera_mat;}
};
struct DirectedLight
{
    float3 dir = float3(0,1,0);
    float intensity = 0;
    float3 color = float3(1,1,1);
    float ambient_q = 1;
    float diffuse_q = 0;
    float specular_q = 0;
    bool has_shadow_map = false;
    float2 shadow_map_size = float2(4096, 4096);
};