#include "vectors.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  dvec2 add(const dvec2 &a, const dvec2 &b)
  {
    dvec2 res{a[0] + b[0], a[1] + b[1]};
    return res;
  }
  dvec2 sub(const dvec2 &a, const dvec2 &b)
  {
    dvec2 res{a[0] - b[0], a[1] - b[1]};
    return res;
  }
  dvec2 mul(const dvec2 &a, const dvec2 &b)
  {
    dvec2 res{a[0] * b[0], a[1] * b[1]};
    return res;
  }
  dvec2 div(const dvec2 &a, const dvec2 &b)
  {
    dvec2 res{a[0] / b[0], a[1] / b[1]};
    return res;
  }
  dvec2 mul(dfloat a, const dvec2 &b)
  {
    dvec2 res{a * b[0], a * b[1]};
    return res;
  }
  dfloat dot(const dvec2 &a, const dvec2 &b)
  {
    return a[0] * b[0] + a[1] * b[1];
  }
  dvec2 normalize(const dvec2 &v)
  {
    dfloat len = CppAD::sqrt(v[0] * v[0] + v[1] * v[1]) + 1e-18;
    dvec2 res;
    res[0] = v[0] / len;
    res[1] = v[1] / len;

    return res;
  }
  dfloat len(const dvec2 &v)
  {
    return CppAD::sqrt(v[0] * v[0] + v[1] * v[1]);
  }


  dvec3 add(const dvec3 &a, const dvec3 &b)
  {
    dvec3 res{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    return res;
  }
  dvec3 sub(const dvec3 &a, const dvec3 &b)
  {
    dvec3 res{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return res;
  }
  dvec3 mul(const dvec3 &a, const dvec3 &b)
  {
    dvec3 res{a[0] * b[0], a[1] * b[1], a[2] * b[2]};
    return res;
  }
  dvec3 div(const dvec3 &a, const dvec3 &b)
  {
    dvec3 res{a[0] / b[0], a[1] / b[1], a[2] / b[2]};
    return res;
  }
  dvec3 mul(dfloat a, const dvec3 &b)
  {
    dvec3 res{a * b[0], a * b[1], a * b[2]};
    return res;
  }
  dfloat dot(const dvec3 &a, const dvec3 &b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }
  dvec3 normalize(const dvec3 &v)
  {
    dfloat len = CppAD::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + 1e-18;
    dvec3 res;
    res[0] = v[0] / len;
    res[1] = v[1] / len;
    res[2] = v[2] / len;

    return res;
  }
  dfloat len(const dvec3 &v)
  {
    return CppAD::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
  dvec3 cross(const dvec3 &a, const dvec3 &b)
  {
    dvec3 res;

    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];

    return res;
  }


  dvec4 add(const dvec4 &a, const dvec4 &b)
  {
    dvec4 res{a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};
    return res;
  }
  dvec4 sub(const dvec4 &a, const dvec4 &b)
  {
    dvec4 res{a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};
    return res;
  }
  dvec4 mul(const dvec4 &a, const dvec4 &b)
  {
    dvec4 res{a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3]};
    return res;
  }
  dvec4 div(const dvec4 &a, const dvec4 &b)
  {
    dvec4 res{a[0] / b[0], a[1] / b[1], a[2] / b[2], a[3] / b[3]};
    return res;
  }
  dvec4 mul(dfloat a, const dvec4 &b)
  {
    dvec4 res{a * b[0], a * b[1], a * b[2], a * b[3]};
    return res;
  }
  dfloat dot(const dvec4 &a, const dvec4 &b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
  }
  dvec4 normalize(const dvec4 &v)
  {
    dfloat len = CppAD::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]) + 1e-18;
    dvec4 res;
    res[0] = v[0] / len;
    res[1] = v[1] / len;
    res[2] = v[2] / len;
    res[3] = v[3] / len;

    return res;
  }
  dfloat len(const dvec4 &v)
  {
    return CppAD::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  }
  dvec4 get_dvec4(const dvec3 &xyz, dfloat w)
  {
    dvec4 res;

    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = w;

    return res;
  }

  dmat43 get_mat43(float *data)
  {
    dmat43 mat;

    for (int i = 0; i < 12; i++)
      mat[i] = data[i];
    
    return mat;
  }

  dmat43 get_mat43(const dvec3 &a, const dvec3 &b, const dvec3 &c, const dvec3 &tr)
  {
    dmat43 mat;

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

    return mat;
  }

  dmat43 ident()
  {
    dmat43 mat;

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

    return mat;
  }

  dmat43 translate(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z)
  {
    dmat43 mat;
    
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

    return mul(input_mat, mat);
  }

  dmat43 translate(const dmat43 &mat, const dvec3 &tr)
  {
    return translate(mat, tr[0], tr[1], tr[2]);
  }

  dmat43 scale(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z)
  {
    dmat43 mat;

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

    return mul(input_mat, mat);
  }

  dmat43 scale(const dmat43 &mat, const dvec3 &sc)
  {
    return scale(mat, sc[0], sc[1], sc[2]);
  }

  dmat43 rotate(const dmat43 &input_mat, const dvec3 &axis, dfloat angle)
  {
    dmat43 mat;
    dvec3 u = normalize(axis);
    dfloat sn = CppAD::sin(angle);
    dfloat cs = CppAD::cos(angle);

    mat[0] = cs + u[0] * u[0] * (1 - cs);
    mat[1] = u[0] * u[1] * (1 - cs) + u[2] * sn;
    mat[2] = u[0] * u[2] * (1 - cs) - u[1] * sn;

    mat[3] = u[0] * u[1] * (1 - cs) - u[2] * sn;
    mat[4] = cs + u[1] * u[1] * (1 - cs);
    mat[5] = u[1] * u[2] * (1 - cs) + u[0] * sn;

    mat[6] = u[0] * u[2] * (1 - cs) + u[1] * sn;
    mat[7] = u[1] * u[2] * (1 - cs) - u[0] * sn;
    mat[8] = cs + u[2] * u[2] * (1 - cs);

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;

    return mul(input_mat, mat);
  }

  dmat43 mul(const dmat43 &a, const dmat43 &b)
  {
    dmat43 mat;
    mat[0] = a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
    mat[1] = a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
    mat[2] = a[2] * b[0] + a[5] * b[1] + a[8] * b[2];

    mat[3] = a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
    mat[4] = a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
    mat[5] = a[2] * b[3] + a[5] * b[4] + a[8] * b[5];

    mat[6] = a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
    mat[7] = a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
    mat[8] = a[2] * b[6] + a[5] * b[7] + a[8] * b[8];

    mat[9] = a[0] * b[9] + a[3] * b[10] + a[6] * b[11] + a[9];
    mat[10] = a[1] * b[9] + a[4] * b[10] + a[7] * b[11] + a[10];
    mat[11] = a[2] * b[9] + a[5] * b[10] + a[8] * b[11] + a[11];
    
    return mat;
  }

  dvec3 mulp(const dmat43 &mat, const dvec3 &vec) // mul (vec.x, vec.y, vec.z, 1)
  {
    dvec3 res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2] + mat[9];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2] + mat[10];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2] + mat[11];
    return res;
  }

  void mulp(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z) // mul (vec.x, vec.y, vec.z, 1)
  {
    dfloat x1 = mat[0] * x + mat[3] * y + mat[6] * z + mat[9];
    dfloat y1 = mat[1] * x + mat[4] * y + mat[7] * z + mat[10];
    dfloat z1 = mat[2] * x + mat[5] * y + mat[8] * z + mat[11];

    x = x1;
    y = y1;
    z = z1;
  }

  dvec3 mulv(const dmat43 &mat, const dvec3 &vec) // mul (vec.x, vec.y, vec.z, 0)
  {
    dvec3 res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2];
    return res;
  }

  void mulv(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z) // mul (vec.x, vec.y, vec.z, 0)
  {
    dfloat x1 = mat[0] * x + mat[3] * y + mat[6] * z;
    dfloat y1 = mat[1] * x + mat[4] * y + mat[7] * z;
    dfloat z1 = mat[2] * x + mat[5] * y + mat[8] * z;

    x = x1;
    y = y1;
    z = z1;
  }

  dvec3 mul4(const dmat43 &mat, const dvec4 &vec)
  {
    dvec3 res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2] + mat[9] * vec[3];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2] + mat[10] * vec[3];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2] + mat[11] * vec[3];
    return res;
  }

  dmat43 transpose3x3(const dmat43 &b)
  {
    dmat43 mat;
    mat[0] = b[0];
    mat[1] = b[3];
    mat[2] = b[6];

    mat[3] = b[1];
    mat[4] = b[4];
    mat[5] = b[7];

    mat[6] = b[2];
    mat[7] = b[5];
    mat[8] = b[8];

    mat[9] = b[9];
    mat[10] = b[10];
    mat[11] = b[11];
    return mat;
  }
  
  dmat43 transposedInverse3x3(const dmat43 &m)
  {
    #define A(i,j) m[3*i + j] 
    dmat43 result;
    dfloat determinant =+A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
                        -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
                        +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    dfloat invdet = 1/(determinant+1e-19);
    result[0] =  (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
    result[3] = -(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
    result[6] =  (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
    result[1] = -(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
    result[4] =  (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
    result[7] = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
    result[2] =  (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
    result[5] = -(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
    result[8] =  (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;
    result[9] =  m[9];
    result[10] = m[10];
    result[11] = m[11];

    return result;
  }

  dmat43 inverse3x4(const dmat43 &m)
  {
    #define A(i,j) m[3*i + j] 
    dmat43 result;
    dfloat determinant =+A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
                        -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
                        +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    dfloat invdet = 1/(determinant+1e-19);
    result[0] =  (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
    result[1] = -(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
    result[2] =  (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
    result[3] = -(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
    result[4] =  (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
    result[5] = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
    result[6] =  (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
    result[7] = -(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
    result[8] =  (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;
    dvec3 inv_tr;
    dvec3 tr{m[9],m[10],m[11]};
    inv_tr = mulv(result, tr);
    result[9] =  -inv_tr[0];
    result[10] = -inv_tr[1];
    result[11] = -inv_tr[2];

    return result;
  }
}