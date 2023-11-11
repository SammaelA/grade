#include "enzyme_test.h"
#include "diff_geometry_generation.h"
#include <cppad/cppad.hpp>
#include <third_party/span.h>
#include <functional>
using dgen::dfloat;

#ifdef USE_CUSTOM_DIFF_RENDER
template <typename my_float>
my_float f4(my_float x, my_float y)
{
  return x+y*y;
}
template <typename my_float>
void f1(my_float *x, my_float *out)
{
  out[0] = 1 + f4(x[0],x[1]);
}

template <typename my_float>
void f2(my_float *x, my_float *out)
{
  dgen::g_vec3<my_float> v1(x[0], x[1], x[2]);
  dgen::g_vec3<my_float> v2(x[3], x[4], x[5]);
  dgen::g_vec3<my_float> res = dgen::cross(v1, v2);
  out[0] = res[0];
  out[1] = res[1];
  out[2] = res[2];
}

template <typename my_float>
void f3(my_float x[1], my_float out[1])
{
  out[0] = 0;
  while (x[0] > 1)
  {
    x[0] -= 1;
    out[0] += 1;
  }
}

void __enzyme_fwddiff(void*, ...);
void __enzyme_autodiff(void*, ...);

#define ENZYME_F(F, x_size, out_size, x, jac) \
{ \
    float out[1024] = {0};\
    float v[1024] = {0};\
    for(int i = 0; i < x_size-1; i++)\
    {\
      v[i] = 1;\
      __enzyme_fwddiff((void*)(F), x, v, enzyme_dupnoneed, out, jac.data() + i*out_size);\
      v[i] = 0;\
    }\
    v[x_size - 1] = 1;\
    __enzyme_fwddiff((void*)(F), x, v, enzyme_dup, out, jac.data() + (x_size - 1)*out_size);\
    logerr("out = %f %f %f", out[0], out[1], out[2]);\
}

#define CPPAD_F(F, x_size, out_size, x, jac) \
  {\
    std::vector<dfloat> X(x_size);\
    std::vector<dfloat> Y(out_size);\
    for (int i=0;i<x_size;i++)\
      X[i] = x[i];\
    CppAD::Independent(X);\
    F<dfloat>(X.data(), Y.data());\
    CppAD::ADFun<float> f(X, Y);\
    std::vector<float> jac_v(out_size * x_size);\
    std::vector<float> X0(x_size);\
    for (int i=0;i<x_size;i++)\
      X0[i] = x[i];\
    jac_v = f.Jacobian(X0);\
    for (int i=0; i<x_size; i++)\
      for (int j=0; j<out_size; j++)\
        jac[i*out_size + j] = jac_v[j*x_size + i];\
  }

#define TEST(test_name, F, x_size, out_size, x) \
{\
  std::vector<float> jac1(x_size*out_size,0);\
  std::vector<float> jac2(x_size*out_size,0);\
  for (int i=0; i<x_size; i++)\
    for (int j=0; j<out_size; j++)\
      jac1[i*out_size + j] = 0;\
  ENZYME_F(F<float>, x_size, out_size, x.data(), jac1);\
  CPPAD_F(F, x_size, out_size, x.data(), jac2);\
  double diff = 0.0;\
  for (int i=0; i<x_size; i++)\
    for (int j=0; j<out_size; j++)\
      diff += abs(jac2[i*out_size + j] - jac1[i*out_size + j]);\
  debug("%s: Jac diff %.2f*1e-9\n",test_name, 1e9*(float)diff);\
  if (diff > 1e-9)\
  {\
  debug("Jacobian enzyme %dx%d\n ---------------- \n",x_size, out_size);\
  for (int i=0; i<x_size; i++)\
  {\
    for (int j=0; j<out_size; j++)\
      debug("%.4f ", jac1[i*out_size + j]);\
    debug("\n");\
  }\
  debug("Jacobian CppAD %dx%d\n ---------------- \n",x_size, out_size);\
  for (int i=0; i<x_size; i++)\
  {\
    for (int j=0; j<out_size; j++)\
      debug("%.4f ", jac2[i*out_size + j]);\
    debug("\n");\
  }\
  }\
}

extern int enzyme_dupnoneed;
extern int enzyme_dup;

void test_enzyme_cppAD()
{
  std::vector<float> x = {1,2};
  TEST("TEST 1", f1, 2, 1, x);

  std::vector<float> v1v2 = {1,0,0,0,1,0};
  TEST("TEST 2", f2, 6, 3, v1v2);

  //enzyme fails to compile here
  //std::vector<float> v3 = {7.5};
  //TEST("TEST 3", f3, 1, 1, v3);
}
#else
void test_enzyme_cppAD()
{
}
#endif