#pragma once
#include "vectors.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  typedef CppAD::AD<float> dfloat;
  typedef g_vec2<dfloat> dvec2;
  typedef g_vec3<dfloat> dvec3;
  typedef g_vec4<dfloat> dvec4;
  typedef g_mat43<dfloat> dmat43;//4 vec3 

  inline dfloat sqrt(dfloat x) { return CppAD::sqrt(x); }
  inline dfloat sin(dfloat x)  { return CppAD::sin(x);  }
  inline dfloat cos(dfloat x)  { return CppAD::cos(x);  }
}
