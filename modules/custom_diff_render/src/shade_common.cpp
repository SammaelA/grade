#include "shade_common.h"
namespace diff_render
{
::std::vector<float> sample_bilinear_clamp(float2 tc, const CPUTexture &tex)
{
  tc = clamp(tc, float2(0,0), float2(1,1));
  tc *= float2(tex.w, tex.h);
  int2 tc0 = clamp(int2(tc), int2(0, 0), int2(tex.w - 1, tex.h - 1));
  int2 tc1 = clamp(int2(tc) + int2(1, 1), int2(0, 0), int2(tex.w - 1, tex.h - 1));
  float2 dtc = tc - float2(tc0);
  const float *p00 = tex.get(tc0.x, tc0.y);
  const float *p01 = tex.get(tc0.x, tc1.y);
  const float *p10 = tex.get(tc1.x, tc0.y);
  const float *p11 = tex.get(tc1.x, tc1.y);

  ::std::vector<float> res(tex.channels, 0);
  for (int i = 0; i < tex.channels; i++)
  {
    res[i] = (1 - dtc.x) * (1 - dtc.y) * p00[i] + (1 - dtc.x) * dtc.y * p01[i] + dtc.x * (1 - dtc.y) * p10[i] + dtc.x * dtc.y * p11[i];
  }

  return res;
}
}