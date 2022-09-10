#pragma once
#include <vector>
#include <array>
namespace CppAD
{
  template <class Base>
  class AD;
}
namespace dgen
{
  typedef CppAD::AD<float> dfloat;
  typedef std::array<dfloat, 2> dvec2;
  typedef std::array<dfloat, 3> dvec3;
  typedef std::array<dfloat, 4> dvec4;
  typedef std::array<dfloat, 12> dmat43;//4 vec3 

  //vec2
  dvec2 add(const dvec2 &a, const dvec2 &b);
  dvec2 sub(const dvec2 &a, const dvec2 &b);
  dvec2 mul(const dvec2 &a, const dvec2 &b);
  dvec2 div(const dvec2 &a, const dvec2 &b);
  dvec2 mul(dfloat a, const dvec2 &b);
  dfloat dot(const dvec2 &a, const dvec2 &b);
  dvec2 normalize(const dvec2 &v);
  dfloat len(const dvec2 &v);

  //vec3
  dvec3 add(const dvec3 &a, const dvec3 &b);
  dvec3 sub(const dvec3 &a, const dvec3 &b);
  dvec3 mul(const dvec3 &a, const dvec3 &b);
  dvec3 div(const dvec3 &a, const dvec3 &b);
  dvec3 mul(dfloat a, const dvec3 &b);
  dfloat dot(const dvec3 &a, const dvec3 &b);
  dvec3 normalize(const dvec3 &v);
  dfloat len(const dvec3 &v);
  dvec3 cross(const dvec3 &a, const dvec3 &b);

  //vec4
  dvec4 add(const dvec4 &a, const dvec4 &b);
  dvec4 sub(const dvec4 &a, const dvec4 &b);
  dvec4 mul(const dvec4 &a, const dvec4 &b);
  dvec4 div(const dvec4 &a, const dvec4 &b);
  dvec4 mul(dfloat a, const dvec4 &b);
  dfloat dot(const dvec4 &a, const dvec4 &b);
  dvec4 normalize(const dvec4 &v);
  dfloat len(const dvec4 &v);
  dvec4 get_dvec4(const dvec3 &xyz, dfloat w);

  //matrices
  dmat43 get_mat43(float *data);
  dmat43 get_mat43(const dvec3 &a, const dvec3 &b, const dvec3 &c, const dvec3 &tr);
  dmat43 ident();
  dmat43 translate(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z);
  dmat43 translate(const dmat43 &mat, const dvec3 &tr);
  dmat43 scale(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z);
  dmat43 scale(const dmat43 &mat, const dvec3 &sc);
  dmat43 rotate(const dmat43 &input_mat, const dvec3 &axis, dfloat angle);
  dmat43 mul(const dmat43 &a, const dmat43 &b);
  dvec3 mulp(const dmat43 &mat, const dvec3 &vec);
  void mulp(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z);
  dvec3 mulv(const dmat43 &mat, const dvec3 &vec);
  void mulv(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z);
  dvec3 mul4(const dmat43 &mat, const dvec4 &vec);
  dmat43 transpose3x3(const dmat43 &b);
  dmat43 transposedInverse3x3(const dmat43 &m);
  dmat43 inverse3x4(const dmat43 &m);
  

}