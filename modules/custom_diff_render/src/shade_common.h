#pragma once
#include "dmodels.h"
namespace diff_render
{
::std::vector<float> sample_bilinear_clamp(float2 tc, const CPUTexture &tex);
}