#include "diff_geometry_generation.h"
#include "common_utils/utility.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  typedef CppAD::AD<float> dfloat;

  typedef dfloat dvec2[2];
  typedef dfloat dvec3[3];
  typedef dfloat dvec4[4];
  typedef dfloat dmat44[16];

  #define FLOAT_PER_VERTEX (3+3+2) //vec3 pos, vec3 norm, vec2 tc
  typedef dfloat dvertex[FLOAT_PER_VERTEX];
  #define MODEL_ARG std::vector<dfloat> &vert
  #define MODEL vert

  inline void get_dvec2(dvec2 res, dfloat x, dfloat y)
  {
    res[0] = x;
    res[1] = y;
  }

  #define __DVEC3_CONSTRUCTOR(t1,t2,t3) \
  inline void get_dvec3(dvec3 res, t1 x, t2 y, t3 z) \
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

  inline void get_dvec4(dvec4 res, dfloat x, dfloat y, dfloat z, dfloat w)
  {
    res[0] = x;
    res[1] = y;
    res[2] = z;
    res[3] = w;
  }

  inline void get_dvec4(dvec4 res, dvec3 xyz, dfloat w)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = w;
  }

  inline void get_dvec4_1(dvec4 res, dvec3 xyz)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = 1;
  }

  inline void get_dvec4_0(dvec4 res, dvec3 xyz)
  {
    res[0] = xyz[0];
    res[1] = xyz[1];
    res[2] = xyz[2];
    res[3] = 0;
  }

  inline void get_dvec4_xyz(dvec4 res, dfloat x, dfloat y, dfloat z)
  {
    res[0] = x;
    res[1] = y;
    res[2] = z;
    res[3] = 1;
  }

  inline void add_dvec3(dvec4 res, const dvec4 a, const dvec4 b)
  {
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];
  }

  inline void add_vertex(MODEL_ARG, int n, const dvec3 pos, const dvec3 norm, const dvec2 tc)
  {
    //int sz = vert.size();
    int sz = n*FLOAT_PER_VERTEX;
    vert.resize(sz + FLOAT_PER_VERTEX);
    vert[sz+0] = pos[0];
    vert[sz+1] = pos[1];
    vert[sz+2] = pos[2];

    vert[sz+3] = norm[0];
    vert[sz+4] = norm[1];
    vert[sz+5] = norm[2]; 

    vert[sz+6] = tc[0];
    vert[sz+7] = tc[1];
  }
/*
  inline void add_triangle(MODEL_ARG, int v1, int v2, int v3)
  {
    int sz = ind.size();
    ind.resize(sz+3);
    ind[sz] = v1;
    ind[sz+1] = v2;
    ind[sz+2] = v3;
  }
*/
  inline int verts(MODEL_ARG)
  {
    return vert.size()/FLOAT_PER_VERTEX;
  }
/*
  inline int triangles(MODEL_ARG)
  {
    return ind.size()/3;
  }
*/
  void test_triangle(MODEL_ARG, std::vector<dfloat> &params)
  {
    dvec3 pos0, pos1, pos2, norm;
    dvec2 tc0, tc1, tc2;
    get_dvec3(pos0, params[0], params[1], params[2]);
    get_dvec3(pos1, params[0]+ params[3], params[1], params[2]);
    get_dvec3(pos2, params[0], params[1], params[2]+ params[4]);
    get_dvec3(norm, 0, 1, 0);
    tc0[0] = 0;
    tc0[1] = 0;
    tc1[0] = 0;
    tc1[1] = 1;
    tc2[0] = 1;
    tc2[1] = 0;

    add_vertex(MODEL, 0, pos0, norm, tc0);
    add_vertex(MODEL, 1, pos1, norm, tc1);
    add_vertex(MODEL, 2, pos2, norm, tc2);
    logerr("triangle");
  }

  void test_cube(MODEL_ARG, std::vector<dfloat> &params)
  {
    
  }

  void Cube(dfloat x0, dfloat y0, dfloat z0,
            dfloat x1, dfloat y1, dfloat z1,
            std::vector<dfloat> &res)
  {
    #define ADD_VEC(x, y, z, i) \
      res[3 * i] = x;           \
      res[3 * i + 1] = y;       \
      res[3 * i + 2] = z;

    ADD_VEC(x0, y0, z0, 0);
    ADD_VEC(x1, y0, z0, 1);
    ADD_VEC(x0, y1, z0, 2);
    ADD_VEC(x1, y1, z0, 3);
    ADD_VEC(x0, y0, z1, 4);
    ADD_VEC(x1, y0, z1, 5);
    ADD_VEC(x0, y1, z1, 6);
    ADD_VEC(x1, y1, z1, 7);
    logerr("cube");
  }
  
  void print_model(const std::vector<float> &res)
  {
    debug("Model\n");
    for (int i=0;i<res.size()/FLOAT_PER_VERTEX;i++)
    {
      int st = i*FLOAT_PER_VERTEX;
      debug("(%.2f %.2f %.2f) (%.2f %.2f %.2f) (%.2f %.2f)\n",
            res[st], res[st+1], res[st+2],
            res[st+3], res[st+4], res[st+5],
            res[st+6], res[st+7]);
    }
  }

  void print_jackobian(const std::vector<float> &jac, int x_n, int y_n)
  {
    std::string names[FLOAT_PER_VERTEX] = {" pos_x", " pos_y", " pos_z", 
                                           "norm_x", "norm_y", "norm_z",
                                           "  tc_x", "  tc_y"};
    debug("Jacobian    ");
    for (int j = 0; j < x_n; j++)
      debug("x_%.2d ", j);
    debugnl();
    for (int i = 0; i < y_n; i++)
    {
      if (i % FLOAT_PER_VERTEX == 0)
      {
        for (int k = 0; k < 13 + 5*x_n; k++)
          debug("-");
        debugnl();
      }
      debug("%s_%.4d ", names[i % FLOAT_PER_VERTEX].c_str(), i);
      for (int j = 0; j < x_n; j++)
        debug("%.2f ", jac[i * x_n + j]);
      debugnl();
    }
  }

  void dgen_test(std::vector<float> &model)
  {
    size_t x_n = 5;
    std::vector<dfloat> X(x_n);
    std::vector<int> inds;
    std::vector<dfloat> Y;

    // declare independent variables and start recording operation sequence
    logerr("gen");
    CppAD::Independent(X);

    //Cube(X[0], X[1], X[2], X[3], X[4], X[5], Y);
    test_triangle(Y, X);
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording
    logerr("gen_finish");
    // compute derivative using operation sequence stored in f
    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    X0[0] = 100;
    X0[1] = 100;
    X0[2] = 0;
    X0[3] = 50;
    X0[4] = 20;

    jac = f.Jacobian(X0); // Jacobian for operation sequence
    res = f.Forward(0, X0);

    print_model(res);
    print_jackobian(jac, x_n, y_n);

    model = res;
  }
}