#include "vectors.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  #define __DVEC2_CONSTRUCTOR(t1,t2) \
  void get_dvec2(dvec2 res, t1 x, t2 y) \
  { \
    res[0] = x;\
    res[1] = y;\
  }
  __DVEC2_CONSTRUCTOR( float,  float)
  __DVEC2_CONSTRUCTOR( float, dfloat)
  __DVEC2_CONSTRUCTOR(dfloat,  float)
  __DVEC2_CONSTRUCTOR(dfloat, dfloat)

  void copy2(dvec2 res, const dvec2 a)
  {
    get_dvec2(res, a[0], a[1]);
  }

  void add2(const dvec2 a, const dvec2 b, dvec2 res)
  {
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
  }

  void sub2(const dvec2 a, const dvec2 b, dvec2 res)
  {
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
  }

  void mul2(const dvec2 a, const dvec2 b, dvec2 res)
  {
    res[0] = a[0] * b[0];
    res[1] = a[1] * b[1];
  }
  void mul2(const dvec2 a, dfloat m, dvec2 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
  }
  void mul2(const dvec2 a, float m, dvec2 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
  }

  void div2(const dvec2 a, const dvec2 b, dvec2 res)
  {
    res[0] = a[0] / b[0];
    res[1] = a[1] / b[1];
  }
  
  dfloat dot2(const dvec2 a, const dvec2 b)
  {
    return a[0]*b[0] + a[1]*b[1];
  }

  void normalize2(dvec2 v)
  {
    dfloat len = CppAD::sqrt(v[0]*v[0] + v[1]*v[1]);
    v[0] /= len;
    v[1] /= len;
  }

  dfloat len2(dvec2 v)
  {
    return CppAD::sqrt(v[0]*v[0] + v[1]*v[1]);
  }

  #define __DVEC3_CONSTRUCTOR(t1,t2,t3) \
  void get_dvec3(dvec3 res, t1 x, t2 y, t3 z) \
  { \
    res[0] = x;\
    res[1] = y;\
    res[2] = z;\
  }

  __DVEC3_CONSTRUCTOR( float,  float,  float)
  __DVEC3_CONSTRUCTOR( float,  float, dfloat)
  __DVEC3_CONSTRUCTOR( float, dfloat,  float)
  __DVEC3_CONSTRUCTOR( float, dfloat, dfloat)
  __DVEC3_CONSTRUCTOR(dfloat,  float,  float)
  __DVEC3_CONSTRUCTOR(dfloat,  float, dfloat)
  __DVEC3_CONSTRUCTOR(dfloat, dfloat,  float)
  __DVEC3_CONSTRUCTOR(dfloat, dfloat, dfloat)

  void copy3(dvec3 res, const dvec3 a)
  {
    get_dvec3(res, a[0], a[1], a[2]);
  }

  void add3(const dvec3 a, const dvec3 b, dvec3 res)
  {
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];
  }

  void sub3(const dvec3 a, const dvec3 b, dvec3 res)
  {
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
    res[2] = a[2] - b[2];
  }

  void mul3(const dvec3 a, const dvec3 b, dvec3 res)
  {
    res[0] = a[0] * b[0];
    res[1] = a[1] * b[1];
    res[2] = a[2] * b[2];
  }
  void mul3(const dvec3 a, dfloat m, dvec3 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
    res[2] = a[2] * m;
  }
  void mul3(const dvec3 a, float m, dvec3 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
    res[2] = a[2] * m;
  }

  void div3(const dvec3 a, const dvec3 b, dvec3 res)
  {
    res[0] = a[0] / b[0];
    res[1] = a[1] / b[1];
    res[2] = a[2] / b[2];
  }
  
  dfloat dot3(const dvec3 a, const dvec3 b)
  {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  void normalize3(dvec3 v)
  {
    dfloat len = CppAD::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
  }

  dfloat len3(dvec3 v)
  {
    return CppAD::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  }

  void cross3(const dvec3 a, const dvec3 b, dvec3 res)
  {
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
  }


  #define __DVEC4_CONSTRUCTOR(t1,t2,t3, t4) \
  void get_dvec4(dvec3 res, t1 x, t2 y, t3 z, t4 w) \
  { \
    res[0] = x;\
    res[1] = y;\
    res[2] = z;\
    res[3] = w;\
  }

  __DVEC4_CONSTRUCTOR( float,   float,  float,  float)
  __DVEC4_CONSTRUCTOR( float,   float,  float, dfloat)
  __DVEC4_CONSTRUCTOR( float,   float, dfloat,  float)
  __DVEC4_CONSTRUCTOR( float,   float, dfloat, dfloat)
  __DVEC4_CONSTRUCTOR( float,  dfloat,  float,  float)
  __DVEC4_CONSTRUCTOR( float,  dfloat,  float, dfloat)
  __DVEC4_CONSTRUCTOR( float,  dfloat, dfloat,  float)
  __DVEC4_CONSTRUCTOR( float,  dfloat, dfloat, dfloat)
  __DVEC4_CONSTRUCTOR(dfloat,   float,  float,  float)
  __DVEC4_CONSTRUCTOR(dfloat,   float,  float, dfloat)
  __DVEC4_CONSTRUCTOR(dfloat,   float, dfloat,  float)
  __DVEC4_CONSTRUCTOR(dfloat,   float, dfloat, dfloat)
  __DVEC4_CONSTRUCTOR(dfloat,  dfloat,  float,  float)
  __DVEC4_CONSTRUCTOR(dfloat,  dfloat,  float, dfloat)
  __DVEC4_CONSTRUCTOR(dfloat,  dfloat, dfloat,  float)
  __DVEC4_CONSTRUCTOR(dfloat,  dfloat, dfloat, dfloat)

  void copy4(dvec4 res, const dvec4 a)
  {
    get_dvec4(res, a[0], a[1], a[2], a[3]);
  }

  void add4(const dvec4 a, const dvec4 b, dvec4 res)
  {
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];
    res[3] = a[3] + b[3];
  }

  void sub4(const dvec4 a, const dvec4 b, dvec4 res)
  {
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
    res[2] = a[2] - b[2];
    res[3] = a[3] - b[3];
  }

  void mul4(const dvec4 a, const dvec4 b, dvec4 res)
  {
    res[0] = a[0] * b[0];
    res[1] = a[1] * b[1];
    res[2] = a[2] * b[2];
    res[3] = a[3] * b[3];
  }
  void mul4(const dvec4 a, dfloat m, dvec4 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
    res[2] = a[2] * m;
    res[3] = a[3] * m;
  }
  void mul4(const dvec4 a, float m, dvec4 res)
  {
    res[0] = a[0] * m;
    res[1] = a[1] * m;
    res[2] = a[2] * m;
    res[3] = a[3] * m;
  }

  void div4(const dvec4 a, const dvec4 b, dvec4 res)
  {
    res[0] = a[0] / b[0];
    res[1] = a[1] / b[1];
    res[2] = a[2] / b[2];
    res[3] = a[3] / b[3];
  }
  
  dfloat dot4(const dvec4 a, const dvec4 b)
  {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
  }

  void normalize4(dvec4 v)
  {
    dfloat len = CppAD::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
    v[3] /= len;
  }

  dfloat len4(dvec4 v)
  {
    return CppAD::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  }

  void get_dvec4(dvec4 res, dvec3 xyz, dfloat w)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = w;
  }

  void get_dvec4_1(dvec4 res, dvec3 xyz)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = 1;
  }

  void get_dvec4_0(dvec4 res, dvec3 xyz)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = 0;
  }

  void get_mat43(dmat43 mat, float *data)
  {
    for (int i=0;i<12;i++)
      mat[i] = data[i];
  }

  void get_mat43(dmat43 mat,const dvec3 a,const dvec3 b, const dvec3 c, const dvec3 tr)
  {
    mat[0] = a[0];
    mat[1] = a[1];
    mat[2] = a[2];

    mat[3] = b[0];
    mat[4] = b[1];
    mat[5] = b[2];

    mat[6] = c[0];
    mat[7] = c[1];
    mat[8] = c[2];

    mat[9] = tr[0];
    mat[10] = tr[1];
    mat[11] = tr[2];
  }

  void ident(dmat43 mat)
  {
    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = 1;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = 1;

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;
  }

  void translate(dmat43 mat, const dfloat x, const dfloat y, const dfloat z)
  {
    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = 1;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = 1;

    mat[9] = x;
    mat[10] = y;
    mat[11] = z;
  }
  
  void translate(dmat43 mat, const dvec3 tr)
  {
    translate(mat, tr[0], tr[1], tr[2]);
  }

  void scale(dmat43 mat, const dfloat x, const dfloat y, const dfloat z)
  {
    mat[0] = x;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = y;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = z;

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;
  }

  void scale(dmat43 mat, const dvec3 tr)
  {
    scale(mat, tr[0], tr[1], tr[2]);
  }

  void rotate(dmat43 mat, const dvec3 axis, const dfloat angle)
  {
    dvec3 u;
    copy3(u, axis);
    normalize3(u);
    dfloat sn = CppAD::sin(angle);
    dfloat cs = CppAD::cos(angle);

    mat[0] = cs + u[0]*u[0]*(1-cs);
    mat[1] = u[0]*u[1]*(1-cs) + u[2]*sn;
    mat[2] = u[0]*u[2]*(1-cs) - u[1]*sn;

    mat[3] = u[0]*u[1]*(1-cs) - u[2]*sn;
    mat[4] = cs + u[1]*u[1]*(1-cs);
    mat[5] = u[1]*u[2]*(1-cs) + u[0]*sn;

    mat[6] = u[0]*u[2]*(1-cs) + u[1]*sn;
    mat[7] = u[1]*u[2]*(1-cs) - u[0]*sn;
    mat[8] = cs + u[2]*u[2]*(1-cs);

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;
  }

  void mul(dmat43 mat, const dmat43 a, const dmat43 b)
  {
    mat[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
    mat[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
    mat[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];

    mat[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
    mat[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
    mat[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];

    mat[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
    mat[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
    mat[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];

    mat[9]  = a[0]*b[9] + a[3]*b[10] + a[6]*b[11] + a[9];
    mat[10] = a[1]*b[9] + a[4]*b[10] + a[7]*b[11] + a[10];
    mat[11] = a[2]*b[9] + a[5]*b[10] + a[8]*b[11] + a[11];
  }

  void mulp(dvec3 res, const dmat43 mat, const dvec3 vec)//mul (vec.x, vec.y, vec.z, 1)
  {
    res[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2] + mat[9];
    res[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2] + mat[10];
    res[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2] + mat[11];
  }

  void mulp(const dmat43 mat, dfloat &x, dfloat &y, dfloat &z)//mul (vec.x, vec.y, vec.z, 1)
  {
    dfloat x1 = mat[0]*x + mat[3]*y + mat[6]*z + mat[9];
    dfloat y1 = mat[1]*x + mat[4]*y + mat[7]*z + mat[10];
    dfloat z1 = mat[2]*x + mat[5]*y + mat[8]*z + mat[11];

    x = x1;
    y = y1;
    z = z1;
  }

  void mulv(dvec3 res, const dmat43 mat, const dvec3 vec)//mul (vec.x, vec.y, vec.z, 0)
  {
    res[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
    res[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
    res[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
  }

  void mulv(const dmat43 mat, dfloat &x, dfloat &y, dfloat &z)//mul (vec.x, vec.y, vec.z, 0)
  {
    dfloat x1 = mat[0]*x + mat[3]*y + mat[6]*z;
    dfloat y1 = mat[1]*x + mat[4]*y + mat[7]*z;
    dfloat z1 = mat[2]*x + mat[5]*y + mat[8]*z;

    x = x1;
    y = y1;
    z = z1;
  }

  void mul4(dvec3 res, const dmat43 mat, const dvec4 vec)
  {
    res[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2] + mat[9]*vec[3];
    res[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2] + mat[10]*vec[3];
    res[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2] + mat[11]*vec[3];
  }
}