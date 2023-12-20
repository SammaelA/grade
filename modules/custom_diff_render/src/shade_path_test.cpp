#include "dmodels.h"
#include "shade_common.h"
#include "utils.h"
namespace diff_render
{
static constexpr int MAX_DEPTH = 8;
static constexpr int RR_DEPTH = 3;

static constexpr float Pi = 3.14159265358979323846;
static constexpr float InvPi = 0.31830988618379067154;
static constexpr float Inv2Pi = 0.15915494309189533577;
static constexpr float Inv4Pi = 0.07957747154594766788;
static constexpr float PiOver2 = 1.57079632679489661923;
static constexpr float PiOver4 = 0.78539816339744830961;
static constexpr float Sqrt2 = 1.41421356237309504880;

inline float CosTheta(const float3 &w) { return w.z; }
inline float Cos2Theta(const float3 &w) { return w.z * w.z; }
inline float AbsCosTheta(const float3 &w) { return ::std::abs(w.z); }
inline float Sin2Theta(const float3 &w) { return ::std::max((float)0, (float)1 - Cos2Theta(w)); }
inline float SinTheta(const float3 &w) { return ::std::sqrt(Sin2Theta(w)); }
inline float TanTheta(const float3 &w) { return SinTheta(w) / CosTheta(w); }
inline float Tan2Theta(const float3 &w) { return Sin2Theta(w) / Cos2Theta(w); }

float rf01()
{
  return ((float)rand())/RAND_MAX;
}
float rf_11()
{
  return 2*rf01() - 1;
}

inline float CosPhi(const float3 &w) 
{
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 1 : clamp(w.x / sinTheta, -1.0f, 1.0f);
}

inline float SinPhi(const float3 &w) 
{
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 0 : clamp(w.y / sinTheta, -1.0f, 1.0f);
}

inline float Cos2Phi(const float3 &w) { return CosPhi(w) * CosPhi(w); }
inline float Sin2Phi(const float3 &w) { return SinPhi(w) * SinPhi(w); }
inline float CosDPhi(const float3 &wa, const float3 &wb) 
{
    float waxy = wa.x * wa.x + wa.y * wa.y;
    float wbxy = wb.x * wb.x + wb.y * wb.y;
    if (waxy == 0 || wbxy == 0)
        return 1;
    return clamp((wa.x * wb.x + wa.y * wb.y) / ::std::sqrt(waxy * wbxy), -1.0f, 1.0f);
}

inline float3 Reflect(const float3 &wo, const float3 &n) { return -1.0f*wo + 2 * dot(wo, n) * n; }

inline bool Refract(const float3 &wi, const float3 &n, float eta, float3 *wt) 
{
  // Compute $\cos \theta_\roman{t}$ using Snell's law
  float cosThetaI = dot(n, wi);
  float sin2ThetaI = ::std::max(float(0), float(1 - cosThetaI * cosThetaI));
  float sin2ThetaT = eta * eta * sin2ThetaI;

  // Handle total internal reflection for transmission
  if (sin2ThetaT >= 1) 
    return false;
  float cosThetaT = ::std::sqrt(1 - sin2ThetaT);
  *wt = eta * (-1.0f*wi) + (eta * cosThetaI - cosThetaT) * n;
  return true;
}

inline bool SameHemisphere(const float3 &w, const float3 &wp) { return w.z * wp.z > 0; }
inline float intensity(const float3 &color) { return 0.2126*color.x + 0.7152*color.y + 0.0722*color.z; }

float2 ConcentricSampleDisk(const float2 &u) 
{
    // Map uniform random numbers to $[-1,1]^2$
    float2 uOffset = 2.f * u - float2(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return float2(0, 0);

    // Apply concentric mapping to point
    float theta, r;
    if (::std::abs(uOffset.x) > ::std::abs(uOffset.y)) {
        r = uOffset.x;
        theta = PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
    }
    return r * float2(::std::cos(theta), ::std::sin(theta));
}

inline float3 CosineSampleHemisphere(const float2 &u) 
{
    float2 d = ConcentricSampleDisk(u);
    float z = ::std::sqrt(::std::max((float)0, 1 - d.x * d.x - d.y * d.y));
    return float3(d.x, d.y, z);
}


struct ExtendedSurfaceInfo
{
  float    t;       // dist from origin ray to surface
  unsigned faceId;  // primitrive id
  float    u;       // first triangle baricentric 
  float    v;       // second triangle baricentric

  float3 pos;       //position of  
  float3 wo;        //normalized direction from surface element to ray start
  float3 n;         //surface normal, also e_z of local coordinate system
  float3 n_geom;    //geometric normal
  float3 e_x;       //x vector of local coordinate system, lays on surface
  float3 e_y;       //y vector of local coordinate system, lays on surface
  ::std::vector<float> sampled_texture;  //sampled texture. Size and content depends on material

  inline float3 WorldToLocal(const float3 &_v) const 
  {
    return float3(dot(_v, e_x), dot(_v, e_y), dot(_v, n));
  }

  inline float3 LocalToWorld(const float3 &_v) const 
  {
    return float3(e_x.x * _v.x + e_y.x * _v.y + n.x * _v.z,
                  e_x.y * _v.x + e_y.y * _v.y + n.y * _v.z,
                  e_x.z * _v.x + e_y.z * _v.y + n.z * _v.z);
  }
};

ExtendedSurfaceInfo get_extended_surface_info(const Scene &scene, const SurfaceInfo &surfInfo, const float3 &pos, const float3 &view_dir)
{
  ExtendedSurfaceInfo res;

  const auto A = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 0);
  const auto B = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 1);
  const auto C = scene.get_index(surfInfo.geomId, surfInfo.instId, surfInfo.primId * 3 + 2);
  const float u = surfInfo.u;
  const float v = surfInfo.v;

  float2 tc = scene.get_tc(A) * (1.0f - u - v) + scene.get_tc(B) * v + u * scene.get_tc(C);
  float3 n = scene.get_norm(A) * (1.0f - u - v) + scene.get_norm(B) * v + u * scene.get_norm(C);
  float3 tangent = scene.get_tang(A) * (1.0f - u - v) + scene.get_tang(B) * v + u * scene.get_tang(C);

  res.t = surfInfo.t;
  res.faceId = surfInfo.primId;
  res.u = surfInfo.u;
  res.v = surfInfo.v;
  res.wo = normalize(-1.0f*view_dir);
  res.n = n;
  res.n_geom = normalize(cross(scene.get_pos(B)-scene.get_pos(A), scene.get_pos(C)-scene.get_pos(A)));
  res.pos = pos + 1e-4*res.n_geom;
  res.e_x = normalize(cross(n,tangent));
  res.e_y = normalize(cross(n,res.e_x)); 
  res.sampled_texture = sample_bilinear_clamp(tc, scene.get_tex(surfInfo.geomId, 0));

  return res;
}

float3 bsdf_val_diffuse(const ExtendedSurfaceInfo &sInfo, const float3 &wo, const float3 &wi)
{
  float3 diffuse = float3(sInfo.sampled_texture[0], sInfo.sampled_texture[1], sInfo.sampled_texture[2]);
  return diffuse / LiteMath::M_PI;
}

float bsdf_pdf_diffuse(const ExtendedSurfaceInfo &sInfo, const float3 &wo, const float3 &wi)
{
  return SameHemisphere(wo, wi) ? AbsCosTheta(wi) / LiteMath::M_PI : 1e9;
}

float3 bsdf_sample_f_diffuse(const ExtendedSurfaceInfo &sInfo, const float3 &wo, float3 *wi, const float2 &u, float *pdf) 
{
  // Cosine-sample the hemisphere, flipping the direction if necessary
  *wi = CosineSampleHemisphere(u);
  if (wo.z < 0) wi->z *= -1;
  *pdf = bsdf_pdf_diffuse(sInfo, wo, *wi);
  return bsdf_val_diffuse(sInfo, wo, *wi);
}

static const ::std::vector<PointLight> pl{PointLight({1,0.8,0.5},2.0, {0,0.7,0.7})};
static ::std::vector<float> pl_cdf;

float3 sample_light(IRayTracer *m_pTracer, const PointLight &light, 
                    const float3 &pos, float3 &light_dir, float &light_pdf, bool &visibility)
{
  float l_sz = 0.1;
  float3 light_pos = light.pos + l_sz*float3(rf_11(),rf_11(),rf_11());
  float3 lv = light_pos - pos;
  float ld = dot(lv, lv);
  light_dir = lv/sqrt(ld);
  visibility = m_pTracer->GetNearestHit(pos, light_dir).primId == unsigned(-1);
  light_pdf = 1;

  return light.intensity*light.color/ld;
}

float3 direct_light(IRayTracer *m_pTracer, const ExtendedSurfaceInfo &sInfo)
{  
  //sample light source with uniform distribution
  int nLights = pl.size();
  const PointLight &light = pl[rand() % nLights];
  float samplePdfInv = nLights;

  //estimate direct light from this source
  float3 L_direct{0,0,0};
  float3 light_dir;
  float light_pdf;
  bool visibility;
  float3 Li = sample_light(m_pTracer, light, sInfo.pos, light_dir, light_pdf, visibility);
  //visibility = true;
  //printf("%f %f %f -- %f %f %f\n", sInfo.pos.x, sInfo.pos.y, sInfo.pos.z, Li.x, Li.y, Li.z);
  if (visibility)
  {
    float3 light_dir_local = sInfo.WorldToLocal(light_dir);
    float3 view_dir_local = sInfo.WorldToLocal(sInfo.wo);

    float3 bsdf_val = bsdf_val_diffuse(sInfo, view_dir_local, light_dir_local)*::std::max(0.0f,dot(sInfo.n, light_dir));
    float bsdf_pdf = bsdf_pdf_diffuse(sInfo, view_dir_local, light_dir_local);

    L_direct = Li*bsdf_val/light_pdf;
    //if (L_direct.x > 0)
    //  printf("%f %f %f -- %f %f %f\n", sInfo.pos.x, sInfo.pos.y, sInfo.pos.z, Li.x, bsdf_val.x, L_direct.x);
  }

  return samplePdfInv*L_direct;
}

template <>
float3 shade<SHADING_MODEL::PATH_TEST>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos)
{
  const float3 ambient = float3(0.5,0.5,0.5);
  float3 L{0,0,0}, beta{1,1,1};
  float3 ray_pos = {0,0,0}, ray_dir = {0,0,0};

  for (int bounce = 0; bounce < MAX_DEPTH; bounce++)
  {
    //find ray-scene intersection
    SurfaceInfo surfInfo;
    if (bounce == 0)
      surfInfo = m_pTracer->CastSingleRay(screen_pos.x, screen_pos.y, &ray_pos, &ray_dir);
    else
      surfInfo = m_pTracer->GetNearestHit(ray_pos, ray_dir);
    
    if (surfInfo.primId == unsigned(-1))
    {
      L += beta*ambient;
      break;
    }
    float3 surf_pos = ray_pos + surfInfo.t*ray_dir;
    ExtendedSurfaceInfo sInfo = get_extended_surface_info(scene, surfInfo, surf_pos, ray_dir);
    L += beta*direct_light(m_pTracer, sInfo);
    //if (L.x > 0)
    //  printf("%f %f %f %f %f -- %f %f %f\n",screen_pos.x, screen_pos.y, surf_pos.x, surf_pos.y, surf_pos.z, L.x, L.y, L.z);
    //float3 light_dir_local = sInfo.WorldToLocal(light_dir);
    float3 view_dir_local = sInfo.WorldToLocal(sInfo.wo);
    float3 next_ray_dir_local;
    float pdf;
    float2 rnd = float2(((float)rand())/RAND_MAX, ((float)rand())/RAND_MAX);

    float3 bsdf_val = bsdf_sample_f_diffuse(sInfo, view_dir_local, &next_ray_dir_local, rnd, &pdf);//*::std::max(0.0f,dot(sInfo.n, light_dir));

    float3 next_ray_dir = sInfo.LocalToWorld(next_ray_dir_local);
    beta *= bsdf_val * ::std::max(0.0f,dot(sInfo.n_geom, next_ray_dir)) / pdf; 
    //printf("eee (%f %f %f) %f %f\n",next_ray_dir_local.x, next_ray_dir_local.y, next_ray_dir_local.z, dot(sInfo.n_geom, next_ray_dir), pdf);

    ray_dir = next_ray_dir;
    ray_pos = sInfo.pos;

    if (bounce > RR_DEPTH) 
    {
      float q = ::std::max((float).05, 1 - intensity(beta));
      if (((float)rand())/RAND_MAX < q)
        break;
      beta /= 1 - q;
    }
  }

  return L;
}

template <>
void shade_grad<SHADING_MODEL::PATH_TEST>(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos,
                                          const float3 val, const AuxData aux, DScene &grad)
{
  shade_grad<SHADING_MODEL::TEXTURE_COLOR>(scene, m_pTracer, screen_pos, val, aux, grad);
}
}