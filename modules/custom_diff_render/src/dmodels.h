#pragma once

#include "LiteMath.h"
#include "dmesh.h"
#include "functions.h"
#include <cstdio>
namespace diff_render
{
#define DEBUG_RENDER 0
constexpr static int  G_DEBUG_VERT_ID = 0;

struct AuxData
{
  const CamInfo* pCamInfo = nullptr;
  Img* debugImages  = nullptr;
  int debugImageNum = 0;
};

template<SHADING_MODEL material>
float3 shade(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos);

template<SHADING_MODEL material>
void shade_grad(const Scene &scene, IRayTracer *m_pTracer, const float2 screen_pos, 
                const float3 val, const AuxData aux, DScene& grad);
}