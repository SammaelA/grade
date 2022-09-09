#pragma once
#include <vector>
namespace CppAD
{
  template <class Base>
  class AD;
}
namespace dgen
{
  typedef CppAD::AD<float> dfloat;
  typedef dfloat dvec2[2];
  typedef dfloat dvec3[3];
  typedef dfloat dvec4[4];
  typedef dfloat dmat43[12];//4 vec3 

  #define FLOAT_PER_VERTEX (3+3+2) //vec3 pos, vec3 norm, vec2 tc
  typedef dfloat dvertex[FLOAT_PER_VERTEX];

  //vec2
  void get_dvec2(dvec2 res, dfloat x, dfloat y);
  void get_dvec2(dvec2 res, dfloat x, float  y);
  void get_dvec2(dvec2 res, float  x, dfloat y);
  void get_dvec2(dvec2 res, float  x, float  y);
  void copy2(dvec2 res, const dvec2 a);
  void add2(const dvec2 a, const dvec2 b, dvec2 res);
  void sub2(const dvec2 a, const dvec2 b, dvec2 res);
  void mul2(const dvec2 a, const dvec2 b, dvec2 res);
  void mul2(const dvec2 a, dfloat m, dvec2 res);
  void mul2(const dvec2 a, float m, dvec2 res);
  void div2(const dvec2 a, const dvec2 b, dvec2 res);
  dfloat dot2(const dvec2 a, const dvec2 b);
  void normalize2(dvec2 v);
  dfloat len2(dvec2 v);
  
  //vec3
  void get_dvec3(dvec3 res, dfloat x, dfloat y, dfloat z);
  void get_dvec3(dvec3 res, dfloat x, dfloat y,  float z);
  void get_dvec3(dvec3 res, dfloat x,  float y, dfloat z);
  void get_dvec3(dvec3 res, dfloat x,  float y,  float z);
  void get_dvec3(dvec3 res,  float x, dfloat y, dfloat z);
  void get_dvec3(dvec3 res,  float x, dfloat y,  float z);
  void get_dvec3(dvec3 res,  float x,  float y, dfloat z);
  void get_dvec3(dvec3 res,  float x,  float y,  float z);
  void copy3(dvec3 res, const dvec3 a);
  void add3(const dvec3 a, const dvec3 b, dvec3 res);
  void sub3(const dvec3 a, const dvec3 b, dvec3 res);
  void mul3(const dvec3 a, const dvec3 b, dvec3 res);
  void mul3(const dvec3 a, dfloat m, dvec3 res);
  void mul3(const dvec3 a, float m, dvec3 res);
  void div3(const dvec3 a, const dvec3 b, dvec3 res);
  dfloat dot3(const dvec3 a, const dvec3 b);
  void normalize3(dvec3 v);
  void normalize3(dfloat &x, dfloat &y, dfloat &z);
  dfloat len3(dvec3 v);
  void cross3(const dvec3 a, const dvec3 b, dvec3 res);

  //vec4
  void get_dvec4(dvec3 res, dfloat x, dfloat y, dfloat z, dfloat w);
  void get_dvec4(dvec3 res, dfloat x, dfloat y, dfloat z,  float w);
  void get_dvec4(dvec3 res, dfloat x, dfloat y,  float z, dfloat w);
  void get_dvec4(dvec3 res, dfloat x, dfloat y,  float z,  float w);
  void get_dvec4(dvec3 res, dfloat x,  float y, dfloat z, dfloat w);
  void get_dvec4(dvec3 res, dfloat x,  float y, dfloat z,  float w);
  void get_dvec4(dvec3 res, dfloat x,  float y,  float z, dfloat w);
  void get_dvec4(dvec3 res, dfloat x,  float y,  float z,  float w);
  void get_dvec4(dvec3 res,  float x, dfloat y, dfloat z, dfloat w);
  void get_dvec4(dvec3 res,  float x, dfloat y, dfloat z,  float w);
  void get_dvec4(dvec3 res,  float x, dfloat y,  float z, dfloat w);
  void get_dvec4(dvec3 res,  float x, dfloat y,  float z,  float w);
  void get_dvec4(dvec3 res,  float x,  float y, dfloat z, dfloat w);
  void get_dvec4(dvec3 res,  float x,  float y, dfloat z,  float w);
  void get_dvec4(dvec3 res,  float x,  float y,  float z, dfloat w);
  void get_dvec4(dvec3 res,  float x,  float y,  float z,  float w);
  void copy4(dvec4 res, const dvec4 a);
  void add4(const dvec4 a, const dvec4 b, dvec4 res);
  void sub4(const dvec4 a, const dvec4 b, dvec4 res);
  void mul4(const dvec4 a, const dvec4 b, dvec4 res);
  void mul4(const dvec4 a, dfloat m, dvec4 res);
  void mul4(const dvec4 a, float m, dvec4 res);
  void div4(const dvec4 a, const dvec4 b, dvec4 res);
  dfloat dot4(const dvec4 a, const dvec4 b);
  void normalize4(dvec4 v);
  dfloat len4(dvec4 v);
  void get_dvec4(dvec4 res, dvec3 xyz, dfloat w);
  void get_dvec4_1(dvec4 res, dvec3 xyz);
  void get_dvec4_0(dvec4 res, dvec3 xyz);
  
  //matrices
  void copy_mat(dmat43 res, dmat43 m);
  void get_mat43(dmat43 mat, float *data);
  void get_mat43(dmat43 mat,const dvec3 a,const dvec3 b, const dvec3 c, const dvec3 tr);
  void ident(dmat43 mat);
  void translate(dmat43 mat, const dfloat x, const dfloat y, const dfloat z);
  void translate(dmat43 mat, const dvec3 tr);
  void scale(dmat43 mat, const dfloat x, const dfloat y, const dfloat z);
  void scale(dmat43 mat, const dvec3 tr);
  void rotate(dmat43 mat, const dvec3 axis, const dfloat angle);
  void mul_mat(dmat43 mat, const dmat43 a, const dmat43 b);
  void mulp(dvec3 res, const dmat43 mat, const dvec3 vec);
  void mulp(const dmat43 mat, dfloat &x, dfloat &y, dfloat &z);
  void mulv(dvec3 res, const dmat43 mat, const dvec3 vec);
  void mulv(const dmat43 mat, dfloat &x, dfloat &y, dfloat &z);
  void mul4(dvec3 res, const dmat43 mat, const dvec4 vec);
  void transpose3x3(dmat43 res, dmat43 m);
  void transpose3x3(dmat43 m);
  void transposedInverse3x3(dmat43 res, dmat43 m);
  void transposedInverse3x3(dmat43 m);
  void inverse3x4(dmat43 res, dmat43 m);
  void inverse3x4(dmat43 m);
}