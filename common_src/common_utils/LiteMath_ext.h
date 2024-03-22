#pragma once
#include "LiteMath/LiteMath.h"
#include <glm/glm.hpp>
#include <cassert>
#include <cstdio>

using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::uint2;
using LiteMath::uint3;
using LiteMath::uint4;
using LiteMath::int2;
using LiteMath::int3;
using LiteMath::int4;
using LiteMath::float4x4;
using LiteMath::float3x3;
using LiteMath::cross;
using LiteMath::dot;
using LiteMath::length;
using LiteMath::normalize;
using LiteMath::to_float4;
using LiteMath::to_float3;

namespace LiteMath
{
  static inline float2 to_float2(float3 f3)         
  { 
    return float2(f3.x, f3.y);
  }
  static inline float to_degrees(float radians)
  {
    return radians * 57.295779513082320876798154814105;
  }
  static inline float to_radians(float degrees)
  {
    return degrees * 0.01745329251994329576923690768489;
  }

  static inline float2 operator-(const float2 v) { return float2{-v.x, -v.y}; }
  static inline float3 operator-(const float3 v) { return float3{-v.x, -v.y, -v.z}; }
  static inline float4 operator-(const float4 v) { return float4{-v.x, -v.y, -v.z, -v.w}; }

  static inline int2 operator-(const int2 v) { return int2{-v.x, -v.y}; }
  static inline int3 operator-(const int3 v) { return int3{-v.x, -v.y, -v.z}; }
  static inline int4 operator-(const int4 v) { return int4{-v.x, -v.y, -v.z, -v.w}; }

  static inline int2 operator%(const int2 a, const int2 b) { return int2{a.x % b.x, a.y % b.y}; }
  static inline int3 operator%(const int3 a, const int3 b) { return int3{a.x % b.x, a.y % b.y, a.z % b.z}; }
  static inline int4 operator%(const int4 a, const int4 b) { return int4{a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w}; }

  static inline uint2 operator%(const uint2 a, const uint2 b) { return uint2{a.x % b.x, a.y % b.y}; }
  static inline uint3 operator%(const uint3 a, const uint3 b) { return uint3{a.x % b.x, a.y % b.y, a.z % b.z}; }
  static inline uint4 operator%(const uint4 a, const uint4 b) { return uint4{a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w}; }

  static inline float4x4 to_float4x4(const float4 col0, const float4 col1, const float4 col2, const float4 col3)
  {
    float4x4 m;
    m.set_col(0, col0);
    m.set_col(1, col1);
    m.set_col(2, col2);
    m.set_col(3, col3);
    return m;
  }

  static inline int2 to_int2(const float2 v){ return int2(v.x, v.y); }
  static inline int3 to_int3(const float3 v){ return int3(v.x, v.y, v.z); }
  static inline int4 to_int4(const float4 v){ return int4(v.x, v.y, v.z, v.w); }

  static inline int2 to_int2(const uint2 v){ return int2(v.x, v.y); }
  static inline int3 to_int3(const uint3 v){ return int3(v.x, v.y, v.z); }
  static inline int4 to_int4(const uint4 v){ return int4(v.x, v.y, v.z, v.w); }

  static inline float2 to_float2(const int2 v){ return float2(v.x, v.y); }
  static inline float3 to_float3(const int3 v){ return float3(v.x, v.y, v.z); }
  static inline float4 to_float4(const int4 v){ return float4(v.x, v.y, v.z, v.w); }

  static inline float2 to_float2(const uint2 v){ return float2(v.x, v.y); }
  static inline float3 to_float3(const uint3 v){ return float3(v.x, v.y, v.z); }
  static inline float4 to_float4(const uint4 v){ return float4(v.x, v.y, v.z, v.w); }

  static inline float2 max(const float2 v, float f){ return float2(std::max(v.x,f), std::max(v.y,f)); }
  static inline float3 max(const float3 v, float f){ return float3(std::max(v.x,f), std::max(v.y,f), std::max(v.z,f)); }
  static inline float4 max(const float4 v, float f){ return float4(std::max(v.x,f), std::max(v.y,f), std::max(v.z,f), std::max(v.w,f)); }

  static inline float2 min(const float2 v, float f){ return float2(std::min(v.x,f), std::min(v.y,f)); }
  static inline float3 min(const float3 v, float f){ return float3(std::min(v.x,f), std::min(v.y,f), std::min(v.z,f)); }
  static inline float4 min(const float4 v, float f){ return float4(std::min(v.x,f), std::min(v.y,f), std::min(v.z,f), std::min(v.w,f)); }

  static inline double sign(double v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }
}