#include "dmodels.h"
#include "shade_common.h"

namespace diff_render
{
template <>
float3 shade<SHADING_MODEL::LAMBERT>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos)
{
  //TODO:move to scene
  const float3 light_dir = normalize(float3(-1,-1,0));
  const float3 light_color = float3(1,0.5,0.2);
  const float3 ambient_light_color = 0.1f*float3(1,1,1);
  const float BIAS = 1e-5;

  float3 ray_pos = {0,0,0}, ray_dir = {0,0,0};
  SurfaceInfo surfInfo = m_pTracer->CastSingleRay(screen_pos.x, screen_pos.y, &ray_pos, &ray_dir);
  if (surfInfo.primId == unsigned(-1))
    return float3(0, 0, 0); // BGCOLOR

  const auto A = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 0);
  const auto B = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 1);
  const auto C = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 2);
  const float u = surfInfo.u;
  const float v = surfInfo.v;

  float2 tc = scene.get_tc(A) * (1.0f - u - v) + scene.get_tc(B) * v + u * scene.get_tc(C);
  float3 n = scene.get_norm(A) * (1.0f - u - v) + scene.get_norm(B) * v + u * scene.get_norm(C);
  auto diffuse_v = sample_bilinear_clamp(tc, scene.get_tex(surfInfo.geomId, 0));
  float3 diffuse = float3(diffuse_v[0], diffuse_v[1], diffuse_v[2]);

  float3 surf_pos = ray_pos + (surfInfo.t-BIAS)*ray_dir;
  float shade = m_pTracer->GetNearestHit(surf_pos, -1.0f*light_dir).primId == unsigned(-1) ? 1 : 0;
  return (shade*light_color*::std::max(0.0f,dot(n,-1.0f*light_dir)) + ambient_light_color)*diffuse;
}

template <>
void shade_grad<SHADING_MODEL::LAMBERT>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos,
                                        const float3 val, const AuxData aux, DScene &grad)
{
  shade_grad<SHADING_MODEL::TEXTURE_COLOR>(scene, m_pTracer, screen_pos, val, aux, grad);
}

template <>
float3 shade<SHADING_MODEL::PHONG>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos)
{
  //TODO:move to scene
  const float3 light_dir = normalize(float3(0,-0.25,-1));
  const float3 light_color = float3(1,0.5,0.2);
  const float3 ambient_light_color = float3(1,1,1);
  const float Ka = 0.1;
  const float Kd = 1;
  const float Ks = 1;
  const int spec_pow = 32;
  const float BIAS = 1e-5;

  float3 ray_pos = {0,0,0}, ray_dir = {0,0,0};
  SurfaceInfo surfInfo = m_pTracer->CastSingleRay(screen_pos.x, screen_pos.y, &ray_pos, &ray_dir);
  if (surfInfo.primId == unsigned(-1))
    return float3(0, 0, 0); // BGCOLOR

  const auto A = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 0);
  const auto B = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 1);
  const auto C = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 2);
  const float u = surfInfo.u;
  const float v = surfInfo.v;

  float2 tc = scene.get_tc(A) * (1.0f - u - v) + scene.get_tc(B) * v + u * scene.get_tc(C);
  float3 n = scene.get_norm(A) * (1.0f - u - v) + scene.get_norm(B) * v + u * scene.get_norm(C);
  auto diffuse_v = sample_bilinear_clamp(tc, scene.get_tex(surfInfo.geomId, 0));
  float3 diffuse = float3(diffuse_v[0], diffuse_v[1], diffuse_v[2]);

  float3 surf_pos = ray_pos + (surfInfo.t-BIAS)*ray_dir;
  float shade = m_pTracer->GetNearestHit(surf_pos, -1.0f*light_dir).primId == unsigned(-1) ? 1 : 0;
  float3 view_dir = ray_dir;
  float3 reflect = light_dir - 2.0f*dot(n,light_dir)*n;
  return (shade*light_color*(Kd*::std::max(0.0f,dot(n,-1.0f*light_dir)) + Ks*pow(::std::max(0.0f,dot(n,reflect)),spec_pow)) + ambient_light_color*Ka)*diffuse;
}

template <>
void shade_grad<SHADING_MODEL::PHONG>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos,
                                     const float3 val, const AuxData aux, DScene &grad)
{
  shade_grad<SHADING_MODEL::TEXTURE_COLOR>(scene, m_pTracer, screen_pos, val, aux, grad);
}

float saturate(float x)
{
  return ::std::max(0.0f, ::std::min(1.0f, x));
}

float G1V(float dotNV, float k)
{
	return 1.0f/(dotNV*(1.0f-k)+k);
}

float LightingFuncGGX(float3 N, float3 V, float3 L, float roughness, float F0)
{
	float alpha = roughness*roughness;

	float3 H = normalize(V+L);

	float dotNL = saturate(dot(N,L));
	float dotNV = saturate(dot(N,V));
	float dotNH = saturate(dot(N,H));
	float dotLH = saturate(dot(L,H));

	float F, D, vis;

	// D
	float alphaSqr = alpha*alpha;
	float pi = 3.14159f;
	float denom = dotNH * dotNH *(alphaSqr-1.0) + 1.0f;
	D = alphaSqr/(pi * denom * denom);

	// F
	float dotLH5 = pow(1.0f-dotLH,5);
	F = F0 + (1.0-F0)*(dotLH5);

	// V
	float k = alpha/2.0f;
	vis = G1V(dotNL,k)*G1V(dotNV,k);

	float specular = dotNL * D * F * vis;
	return specular;
}

template <>
float3 shade<SHADING_MODEL::GGX>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos)
{
  //TODO:move to scene
  const float3 light_dir = normalize(float3(0.15,0,-1));
  const float3 light_color = float3(1,0.5,0.2);
  const float3 ambient_light_color = float3(1,1,1);
  const float Ka = 0.01;
  const float BIAS = 1e-5;

  float3 ray_pos = {0,0,0}, ray_dir = {0,0,0};
  SurfaceInfo surfInfo = m_pTracer->CastSingleRay(screen_pos.x, screen_pos.y, &ray_pos, &ray_dir);
  if (surfInfo.primId == unsigned(-1))
    return float3(0, 0, 0); // BGCOLOR

  const auto A = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 0);
  const auto B = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 1);
  const auto C = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 2);
  const float u = surfInfo.u;
  const float v = surfInfo.v;

  float2 tc = scene.get_tc(A) * (1.0f - u - v) + scene.get_tc(B) * v + u * scene.get_tc(C);
  float3 n = scene.get_norm(A) * (1.0f - u - v) + scene.get_norm(B) * v + u * scene.get_norm(C);
  auto diffuse_v = sample_bilinear_clamp(tc, scene.get_tex(surfInfo.geomId, 0));
  float3 diffuse = float3(diffuse_v[0], diffuse_v[1], diffuse_v[2]);

  float3 surf_pos = ray_pos + (surfInfo.t-BIAS)*ray_dir;
  float shade = m_pTracer->GetNearestHit(surf_pos, -1.0f*light_dir).primId == unsigned(-1) ? 1 : 0;
  float3 view_dir = ray_dir;
  float3 reflect = light_dir - 2.0f*dot(n,light_dir)*n;
  float ggx_res = LightingFuncGGX(n, -1.0f*view_dir, -1.0f*light_dir, 0.2, 0.8);
  return (shade*light_color*(::std::max(0.0f,dot(n,-1.0f*light_dir)) + ggx_res) + Ka*ambient_light_color)*diffuse;
}

template <>
void shade_grad<SHADING_MODEL::GGX>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos,
                                    const float3 val, const AuxData aux, DScene &grad)
{
  shade_grad<SHADING_MODEL::TEXTURE_COLOR>(scene, m_pTracer, screen_pos, val, aux, grad);
}
}