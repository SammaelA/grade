#include "drender.h"
namespace diff_render
{
// build a discrete CDF using edge length
Sampler build_edge_sampler(const Scene &scene, const ::std::vector<Edge> &edges) 
{
  ::std::vector<float> pmf;
  ::std::vector<float> cdf;
  pmf.reserve(edges.size());
  cdf.reserve(edges.size() + 1);
  cdf.push_back(0);
  for (auto edge : edges) {
      auto v0 = scene.get_pos(edge.v0);
      auto v1 = scene.get_pos(edge.v1);
      pmf.push_back(length(v1 - v0));
      cdf.push_back(pmf.back() + cdf.back());
  }
  auto length_sum = cdf.back();
  ::std::for_each(pmf.begin(), pmf.end(), [&](float &p) {p /= length_sum;});
  ::std::for_each(cdf.begin(), cdf.end(), [&](float &p) {p /= length_sum;});
  return Sampler{pmf, cdf};
}

// binary search for inverting the CDF in the sampler
int sample(const Sampler &sampler, const float u) 
{
  auto cdf = sampler.cdf;
  return clamp(::std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin() - 1, 0, cdf.size() - 2);
}

inline void edge_grad(const Scene &scene, const Edge &e, const float2 d_v0, const float2 d_v1, const AuxData aux,
                      ::std::vector<::std::vector<GradReal>> &d_pos, ::std::vector<::std::vector<GradReal>> &d_tr)
{
  float3 v0_d[2] = {{0, 0, 0}, {0, 0, 0}};
  float3 v1_d[2] = {{0, 0, 0}, {0, 0, 0}};

  float3 v0_3d = scene.get_pos(e.v0);
  float3 v1_3d = scene.get_pos(e.v1);

  VS_X_grad(v0_3d.M, *(aux.pCamInfo), v0_d[0].M);
  VS_Y_grad(v0_3d.M, *(aux.pCamInfo), v0_d[1].M);
  VS_X_grad(v1_3d.M, *(aux.pCamInfo), v1_d[0].M);
  VS_Y_grad(v1_3d.M, *(aux.pCamInfo), v1_d[1].M);

  const float dv0_dx = v0_d[0].x * d_v0.x; // + v0_dx.y*d_v0.y; ==> 0
  const float dv0_dy = v0_d[1].y * d_v0.y; // + v0_dy.x*d_v0.x; ==> 0
  const float dv0_dz = (v0_d[0].z * d_v0.x + v0_d[1].z * d_v0.y);

  const float dv1_dx = v1_d[0].x * d_v1.x; // + v1_dx.y*d_v1.y; ==> 0
  const float dv1_dy = v1_d[1].y * d_v1.y; // + v1_dy.x*d_v1.x; ==> 0
  const float dv1_dz = (v1_d[0].z * d_v1.x + v1_d[1].z * d_v1.y);

  //v0 = scene.get_index(mesh_id, instance_id, primitive_id)
  //v0_3d = mesh[mesh_id].transforms[instance_id] * mesh[mesh_id].vertices[primitive_id]
  //                       |                                     |
  //                      tr                                     v
  //we should calculate  both  d_v0/d_tr(d_v1/d_tr) and d_v0/d_v(d_v1/d_v)

  int tr_n = e.instance_n;
  float4x4 tr = scene.get_transform_inv(0)[tr_n];
  float3 v0_orig = scene.get_pos_orig(e.v0);
  float3 v1_orig = scene.get_pos_orig(e.v1);
  //logerr("set %d %f %f %f %f %f %f",v0,tr(0,0),dv0_dx,tr(0,1),dv0_dy,tr(0,2),dv0_dz);

  d_pos[e.mesh_n][e.vn0 * 3 + 0] += tr(0,0)*dv0_dx + tr(0,1)*dv0_dy + tr(0,2)*dv0_dz;
  d_pos[e.mesh_n][e.vn0 * 3 + 1] += tr(1,0)*dv0_dx + tr(1,1)*dv0_dy + tr(1,2)*dv0_dz;
  d_pos[e.mesh_n][e.vn0 * 3 + 2] += tr(2,0)*dv0_dx + tr(2,1)*dv0_dy + tr(2,2)*dv0_dz;

  d_pos[e.mesh_n][e.vn1 * 3 + 0] += tr(0,0)*dv1_dx + tr(0,1)*dv1_dy + tr(0,2)*dv1_dz;
  d_pos[e.mesh_n][e.vn1 * 3 + 1] += tr(1,0)*dv1_dx + tr(1,1)*dv1_dy + tr(1,2)*dv1_dz;
  d_pos[e.mesh_n][e.vn1 * 3 + 2] += tr(2,0)*dv1_dx + tr(2,1)*dv1_dy + tr(2,2)*dv1_dz;

  d_tr[e.mesh_n][12*tr_n + 0] += v0_orig.x*dv0_dx + v1_orig.x*dv1_dx;
  d_tr[e.mesh_n][12*tr_n + 1] += v0_orig.y*dv0_dx + v1_orig.y*dv1_dx;
  d_tr[e.mesh_n][12*tr_n + 2] += v0_orig.z*dv0_dx + v1_orig.z*dv1_dx;
  d_tr[e.mesh_n][12*tr_n + 3] +=         1*dv0_dx +         1*dv1_dx;

  d_tr[e.mesh_n][12*tr_n + 4] += v0_orig.x*dv0_dy + v1_orig.x*dv1_dy;
  d_tr[e.mesh_n][12*tr_n + 5] += v0_orig.y*dv0_dy + v1_orig.y*dv1_dy;
  d_tr[e.mesh_n][12*tr_n + 6] += v0_orig.z*dv0_dy + v1_orig.z*dv1_dy;
  d_tr[e.mesh_n][12*tr_n + 7] +=         1*dv0_dy +         1*dv1_dy;

  d_tr[e.mesh_n][12*tr_n + 8] += v0_orig.x*dv0_dz + v1_orig.x*dv1_dz;
  d_tr[e.mesh_n][12*tr_n + 9] += v0_orig.y*dv0_dz + v1_orig.y*dv1_dz;
  d_tr[e.mesh_n][12*tr_n + 10]+= v0_orig.z*dv0_dz + v1_orig.z*dv1_dz;
  d_tr[e.mesh_n][12*tr_n + 11]+=         1*dv0_dz +         1*dv1_dz;
}

::std::shared_ptr<IDiffRender> MakeDifferentialRenderer(const DiffRenderSettings &settings)
{
  switch (settings.mode)
  {
  case SHADING_MODEL::SILHOUETTE:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::SILHOUETTE>>();
    impl->init(settings);
    return impl;
    }
    break;    
  case SHADING_MODEL::VERTEX_COLOR:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::VERTEX_COLOR>>();
    impl->init(settings);
    return impl;
    }
    break;
  case SHADING_MODEL::TEXTURE_COLOR:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::TEXTURE_COLOR>>();
    impl->init(settings);
    return impl;
    }
    break;
  case SHADING_MODEL::LAMBERT:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::LAMBERT>>();
    impl->init(settings);
    return impl;
    }
    break;
  case SHADING_MODEL::PHONG:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::PHONG>>();
    impl->init(settings);
    return impl;
    }
    break;
  case SHADING_MODEL::GGX:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::GGX>>();
    impl->init(settings);
    return impl;
    }
    break;
  case SHADING_MODEL::PATH_TEST:
    {
    auto impl = ::std::make_shared<DiffRender<SHADING_MODEL::PATH_TEST>>();
    impl->init(settings);
    return impl;
    }
    break;
  default:
    assert(false);
    break;
  }
  return nullptr;
}
}