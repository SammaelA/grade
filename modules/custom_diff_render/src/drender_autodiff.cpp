#include "drender.h"

void __enzyme_autodiff(void*, ...);
void __enzyme_fwddiff(void*, ...);
extern int enzyme_dupnoneed;

namespace diff_render
{
inline void translate(float mat[DMesh::TRANSFORM_SIZE], float x, float y, float z) 
{
  mat[0*4 + 3] += x;
  mat[1*4 + 3] += y;
  mat[2*4 + 3] += z;
}

//out = A*B
inline void mul(float a[DMesh::TRANSFORM_SIZE], float b[DMesh::TRANSFORM_SIZE], float out[DMesh::TRANSFORM_SIZE])
{
  /*
  0  1  2  3
  4  5  6  7
  8  9 10 11
  */
  out[0] = a[0]*b[0] + a[1]*b[4] + a[2]*b[8];
  out[1] = a[0]*b[1] + a[1]*b[5] + a[2]*b[9];
  out[2] = a[0]*b[2] + a[1]*b[6] + a[2]*b[10];
  out[3] = a[0]*b[3] + a[1]*b[7] + a[2]*b[11] + a[3];

  out[4] = a[4]*b[0] + a[5]*b[4] + a[6]*b[8];
  out[5] = a[4]*b[1] + a[5]*b[5] + a[6]*b[9];
  out[6] = a[4]*b[2] + a[5]*b[6] + a[6]*b[10];
  out[7] = a[4]*b[3] + a[5]*b[7] + a[6]*b[11] + a[7];

  out[8] = a[8]*b[0] + a[9]*b[4] + a[10]*b[8];
  out[9] = a[8]*b[1] + a[9]*b[5] + a[10]*b[9];
  out[10]= a[8]*b[2] + a[9]*b[6] + a[10]*b[10];
  out[11]= a[8]*b[3] + a[9]*b[7] + a[10]*b[11] + a[11];
}

inline void scale(float s, float out[DMesh::TRANSFORM_SIZE])
{
  out[0] = s   ; out[1] = 0.0f; out[2] = 0.0f; out[3] = 0.0f;
  out[4] = 0.0f; out[5] = s   ; out[6] = 0.0f; out[7] = 0.0f;
  out[8] = 0.0f; out[9] = 0.0f; out[10]= s   ; out[11]= 0.0f;  
}

inline void rotate_x(float r, float out[DMesh::TRANSFORM_SIZE])
{
  float s = sinf(r);
  float c = cosf(r);
  out[0] = 1   ; out[1] = 0.0f; out[2] = 0.0f; out[3] = 0.0f;
  out[4] = 0.0f; out[5] = c   ; out[6] = -s  ; out[7] = 0.0f;
  out[8] = 0.0f; out[9] = s   ; out[10]= c   ; out[11]= 0.0f;  
}

inline void rotate_y(float r, float out[DMesh::TRANSFORM_SIZE])
{
  float s = sinf(r);
  float c = cosf(r);
  out[0] = c   ; out[1] = 0.0f; out[2] = s   ; out[3] = 0.0f;
  out[4] = 0.0f; out[5] = 1   ; out[6] = 0.0f; out[7] = 0.0f;
  out[8] = -s  ; out[9] = 0.0f; out[10]= c   ; out[11]= 0.0f;  
}

inline void rotate_z(float r, float out[DMesh::TRANSFORM_SIZE])
{
  float s = sinf(r);
  float c = cosf(r);
  out[0] = c   ; out[1] = -s  ; out[2] = 0.0f; out[3] = 0.0f;
  out[4] = s   ; out[5] = c   ; out[6] = 0.0f; out[7] = 0.0f;
  out[8] = 0.0f; out[9] = 0.0f; out[10]= 1   ; out[11]= 0.0f;  
}

void dmat_dtransforms_jac_f(float x[DMesh::RESTRICTED_TRANSFORM_SIZE], float out[DMesh::TRANSFORM_SIZE]) 
{
  float sc   [DMesh::TRANSFORM_SIZE];
  float rot_x[DMesh::TRANSFORM_SIZE];
  float rot_y[DMesh::TRANSFORM_SIZE];
  float rot_z[DMesh::TRANSFORM_SIZE];

  scale(x[6], sc);
  rotate_x(x[3], rot_x);
  rotate_y(x[4], rot_y);
  rotate_z(x[5], rot_z);

  mul(rot_x, sc, out);//out = rot_x*scale
  mul(rot_y, out, sc);//sc = rot_y*rot_x*scale
  mul(rot_z, sc, out);//out = rot_z*rot_y*rot_x*scale

  translate(out, x[0], x[1], x[2]);//out = translate*rot_z*rot_y*rot_x*scale
}
}
void dmat_dtransforms_jac(float const mat[diff_render::DMesh::RESTRICTED_TRANSFORM_SIZE], float jac[diff_render::DMesh::RESTRICTED_TRANSFORM_SIZE][diff_render::DMesh::TRANSFORM_SIZE])
{
    float out[diff_render::DMesh::TRANSFORM_SIZE] = {0};
    float v[diff_render::DMesh::RESTRICTED_TRANSFORM_SIZE] = {0};
    for(int i = 0; i < diff_render::DMesh::RESTRICTED_TRANSFORM_SIZE; i++) 
    {
      v[i] = 1;
      __enzyme_fwddiff((void*)diff_render::dmat_dtransforms_jac_f, mat, v, enzyme_dupnoneed, out, jac[i]);
      v[i] = 0;
    }
}