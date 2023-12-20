#include "raytrace.h"
using LiteMath::dot;
using LiteMath::sign;
using LiteMath::cross;
using LiteMath::float4x4;
using LiteMath::float3;
using LiteMath::normalize;
using LiteMath::inverse4x4;
using LiteMath::to_float3;
namespace diff_render
{
//#include <iostream>

//static inline float BarU( const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3])
//{
//  const float edge1X = B[0] - A[0];
//  const float edge1Y = B[1] - A[1];
//  const float edge1Z = B[2] - A[2];
//
//  const float edge2X = C[0] - A[0];
//  const float edge2Y = C[1] - A[1];
//  const float edge2Z = C[2] - A[2];
//  
//  const float pvecZ = ray_dir[0]*edge2Y - ray_dir[1]*edge2X;
//  const float pvecX = ray_dir[1]*edge2Z - ray_dir[2]*edge2Y;
//  const float pvecY = ray_dir[2]*edge2X - ray_dir[0]*edge2Z;
//
//  const float tvecX  = ray_pos[0] - A[0];
//  const float tvecY  = ray_pos[1] - A[1];
//  const float tvecZ  = ray_pos[2] - A[2];
//
//  const float qvecZ  = tvecX*edge1Y - tvecY*edge1X;
//  const float qvecX  = tvecY*edge1Z - tvecZ*edge1Y;
//  const float qvecY  = tvecZ*edge1X - tvecX*edge1Z;
//
//  const float e1dp   = edge1X*pvecX + edge1Y*pvecY + edge1Z*pvecZ;
//  const float signv  = sign(e1dp); // put 1.0 to enable triangle clippin
//
//  return signv*(qvecX*ray_dir[0] + qvecY*ray_dir[1] + qvecZ*ray_dir[2])/(e1dp*signv);
//}

//static inline float BarV(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3])
//{
//  const float edge1X = B[0] - A[0];
//  const float edge1Y = B[1] - A[1];
//  const float edge1Z = B[2] - A[2];
//
//  const float edge2X = C[0] - A[0];
//  const float edge2Y = C[1] - A[1];
//  const float edge2Z = C[2] - A[2];
//
//  const float pvecZ = ray_dir[0]*edge2Y - ray_dir[1]*edge2X;
//  const float pvecX = ray_dir[1]*edge2Z - ray_dir[2]*edge2Y;
//  const float pvecY = ray_dir[2]*edge2X - ray_dir[0]*edge2Z;
//
//  const float tvecX  = ray_pos[0] - A[0];
//  const float tvecY  = ray_pos[1] - A[1];
//  const float tvecZ  = ray_pos[2] - A[2];
//
//  const float qvecZ  = tvecX*edge1Y - tvecY*edge1X;
//  const float qvecX  = tvecY*edge1Z - tvecZ*edge1Y;
//  const float qvecY  = tvecZ*edge1X - tvecX*edge1Z;
//
//  const float e1dp   = edge1X*pvecX + edge1Y*pvecY + edge1Z*pvecZ;
//  const float signv  = sign(e1dp); // put 1.0 to enable triangle clippin
//
//  return signv*(tvecX*pvecX + tvecY*pvecY + tvecZ*pvecZ)/(signv*e1dp); 
//}

struct BruteForce3D : public IRayTracer
{
  BruteForce3D(){}
  ~BruteForce3D() override {}

  void Init(const Scene* pScene) override 
  {
    m_pScene = pScene;
    //::std::cout << "[BruteForce3D]: Init done" << ::std::endl;
    int max_instances = 0;
    for (const auto &t : pScene->get_transforms())
      max_instances = ::std::max(max_instances, (int)t.size());
    ::std::vector<int> instance_ids;
    for (int i=0;i<max_instances;i++)
      instance_ids.push_back(i);
    pScene->set_instance_id_mapping(instance_ids);
  }

  float3 GetCameraPos() const override
  {
    return m_camPos;
  }

  void SetCamera(const CamInfo& cam) override
  {
    m_ProjInv      = inverse4x4(cam.mProj);
    m_worldViewInv = inverse4x4(cam.mWorldView);
    m_fwidth       = cam.width;
    m_fheight      = cam.height;
    m_camPos = float3(cam.mWorldView.get_col(3).x, cam.mWorldView.get_col(3).y, cam.mWorldView.get_col(3).z);
  }

  SurfaceInfo GetNearestHit(float3 rayPos, float3 rayDir, float tNear = 0.0f, float tFar = 1e9f) override
  {
    SurfaceInfo hit;
    hit.geomId = unsigned(-1);
    hit.instId = unsigned(-1);
    hit.primId = unsigned(-1);
    hit.u      = 0.0f;
    hit.v      = 0.0f;
    hit.t      = tFar; // tFar
    
    for (int i=0;i<m_pScene->get_meshes().size();i++)
    {
      for (int j=0;j<m_pScene->get_transform(i).size();j++)
      {
        float4x4 transform = m_pScene->get_transform(i)[j];
        for (size_t triAddress = 0; triAddress < m_pScene->get_mesh(i).indices.size(); triAddress += 3)
        { 
          const uint A = m_pScene->get_index(i, j, triAddress + 0);
          const uint B = m_pScene->get_index(i, j, triAddress + 1);
          const uint C = m_pScene->get_index(i, j, triAddress + 2);
        
          const float3 A_pos = transform*m_pScene->get_pos(A);
          const float3 B_pos = transform*m_pScene->get_pos(B);
          const float3 C_pos = transform*m_pScene->get_pos(C);
        
          const float3 edge1 = B_pos - A_pos;
          const float3 edge2 = C_pos - A_pos;
          const float3 pvec  = cross(rayDir, edge2);
          const float3 tvec  = rayPos - A_pos;
          const float3 qvec  = cross(tvec, edge1);
          const float  e1dpv = dot(edge1, pvec);
          const float  signv = sign(e1dpv);                 // put 1.0 to enable triangle clipping
          const float invDet = signv / ::std::max(signv*e1dpv, 1e-6f);

          const float v = dot(tvec, pvec)*invDet;
          const float u = dot(qvec, rayDir)*invDet;
          const float t = dot(edge2, qvec)*invDet;
        
          if (v > 0.0f && u > 0.0f && (u + v < 1.0f) && t > tNear && t < hit.t)
          {
            //const float u2 = BarU(ray_pos.M, ray_dir.M, A_pos.M, B_pos.M, C_pos.M);
            //const float v2 = BarV(ray_pos.M, ray_dir.M, A_pos.M, B_pos.M, C_pos.M);
            hit.t      = t;
            hit.geomId = i;
            hit.instId = j;
            hit.primId = triAddress/3;
            hit.u      = u;    // v2
            hit.v      = v;    // v1
          }
        }
      }
    }
  
    return hit;
  }

  SurfaceInfo CastSingleRay(float x, float y, float3* outPos, float3* outDir) override
  {
    const float2 screen_pos(x,y);
  
    float3 ray_pos = float3(0,0,0);
    float3 ray_dir = EyeRayDirNormalized(x/m_fwidth, y/m_fheight, m_ProjInv);
    float  tNear   = 0.0f;

    transform_ray3f(m_worldViewInv, &ray_pos, &ray_dir);

    if(outPos != nullptr)
      *outPos = ray_pos;
    if(outDir != nullptr)
      *outDir = ray_dir;

    return GetNearestHit(ray_pos, ray_dir);
  }

  const Scene* m_pScene = nullptr;
  float4x4 m_ProjInv;
  float4x4 m_worldViewInv;
  float3 m_camPos;
  float m_fwidth, m_fheight;

};

#ifdef USE_EMBREE
::std::shared_ptr<IRayTracer> MakeEmbreeRT3D();
::std::shared_ptr<IRayTracer> MakeRayTracer3D(const char* className) { return MakeEmbreeRT3D(); }
#else
::std::shared_ptr<IRayTracer> MakeRayTracer3D(const char* className) { return ::std::make_shared<BruteForce3D>(); }
#endif
}