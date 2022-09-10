#include "diff_geometry_generation.h"
#include "common_utils/utility.h"
#include "graphics_utils/modeling.h"
#include <cppad/cppad.hpp>
#include "tinyEngine/engine.h"
namespace dgen
{
  #define FLOAT_PER_VERTEX (3+3+2) //vec3 pos, vec3 norm, vec2 tc

  inline void add_vertex(std::vector<dfloat> &vert, int n, const dvec3 pos, const dvec3 norm, const dvec2 tc)
  {
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

  inline int verts(std::vector<dfloat> &vert)
  {
    return vert.size()/FLOAT_PER_VERTEX;
  }

  inline std::vector<float> get_triangle()
  {
    std::vector<float> res = {
      0,0,0, 0,0,1, 0,0,
      1,0,0, 0,0,1, 0,1,
      0,1,0, 0,0,1, 1,0
    };
    return res;
  }

  inline std::vector<float> get_cube_expl()
  {
    std::vector<float> res = {
      0.00,0.00,0.00, -1.00,-0.00,-0.00, 0.00,0.00,
      0.00,1.00,0.00, -1.00,-0.00,-0.00, 0.00,1.00,
      0.00,0.00,1.00, -1.00,-0.00,-0.00, 1.00,0.00,
      0.00,1.00,1.00, -1.00,-0.00,-0.00, 1.00,1.00,
      0.00,0.00,1.00, -1.00,-0.00,-0.00, 1.00,0.00,
      0.00,1.00,0.00, -1.00,-0.00,-0.00, 0.00,1.00,

      0.00,1.00,0.00, 0.00,1.00, 0.00, 0.00,1.00,
      0.00,1.00,1.00, 0.00,1.00, 0.00, 1.00,1.00,
      1.00,1.00,0.00, 0.00,1.00, 0.00, 0.00,1.00,
      1.00,1.00,1.00, 0.00,1.00, 0.00, 1.00,1.00,
      1.00,1.00,0.00, 0.00,1.00, 0.00, 0.00,1.00,
      0.00,1.00,1.00, 0.00,1.00, 0.00, 1.00,1.00,

      1.00,1.00,1.00, 1.00,-0.00,-0.00, 1.00,1.00,
      1.00,0.00,1.00, 1.00,-0.00,-0.00, 1.00,0.00,
      1.00,1.00,0.00, 1.00,-0.00,-0.00, 0.00,1.00,
      1.00,0.00,0.00, 1.00,-0.00,-0.00, 0.00,0.00,
      1.00,1.00,0.00, 1.00,-0.00,-0.00, 0.00,1.00,
      1.00,0.00,1.00, 1.00,-0.00,-0.00, 1.00,0.00,
      
      0.00,0.00,0.00, 0.00,-1.00,0.00, 0.00,0.00,
      0.00,0.00,1.00, 0.00,-1.00,0.00, 1.00,0.00,
      1.00,0.00,0.00, 0.00,-1.00,0.00, 0.00,0.00,
      1.00,0.00,1.00, 0.00,-1.00,0.00, 1.00,0.00,
      0.00,0.00,1.00, 0.00,-1.00,0.00, 1.00,0.00,
      1.00,0.00,0.00, 0.00,-1.00,0.00, 0.00,0.00,

      0.00,0.00,0.00, 0.00,-0.00,-1.00, 0.00,0.00,
      1.00,0.00,0.00, 0.00,-0.00,-1.00, 0.00,0.00,
      0.00,1.00,0.00, 0.00,-0.00,-1.00, 0.00,1.00,
      1.00,1.00,0.00, 0.00,-0.00,-1.00, 0.00,1.00,
      0.00,1.00,0.00, 0.00,-0.00,-1.00, 0.00,1.00,
      1.00,0.00,0.00, 0.00,-0.00,-1.00, 0.00,0.00,

      0.00,0.00,1.00, 0.00,-0.00,1.00, 1.00,0.00,
      1.00,0.00,1.00, 0.00,-0.00,1.00, 1.00,0.00,
      0.00,1.00,1.00, 0.00,-0.00,1.00, 1.00,1.00,
      1.00,1.00,1.00, 0.00,-0.00,1.00, 1.00,1.00,
      0.00,1.00,1.00, 0.00,-0.00,1.00, 1.00,1.00,
      1.00,0.00,1.00, 0.00,-0.00,1.00, 1.00,0.00
    };
    return res;
  }

  void model_to_simple_model(Mesh *m, std::vector<float> &s_model)
  {
    s_model.resize(m->indices.size()*FLOAT_PER_VERTEX);
    int pos = 0;
    bool have_normals = (m->normals.size() == m->positions.size());
    for (int ind : m->indices)
    {
      s_model[pos] = m->positions[3*ind];
      s_model[pos+1] = m->positions[3*ind+1];
      s_model[pos+2] = m->positions[3*ind+2];
      
      if (have_normals)
      {
        s_model[pos+3] = m->normals[3*ind];
        s_model[pos+4] = m->normals[3*ind+1];
        s_model[pos+5] = m->normals[3*ind+2];
      }
      else
      {
        s_model[pos+3] = 1;
        s_model[pos+4] = 0;
        s_model[pos+5] = 0;
      }
      if (m->colors.size()/4 >= m->positions.size()/3)
      {
        s_model[pos+6] = m->colors[4*ind];
        s_model[pos+7] = m->colors[4*ind+1];
      }
      else
      {
        s_model[pos+6] = 0;
        s_model[pos+7] = 0;
      }
      pos+=FLOAT_PER_VERTEX;
    }
  }

  std::vector<float> get_cube()
  {
    Box b = Box(glm::vec3(0,0,0), glm::vec3(1,0,0), glm::vec3(0,1,0), glm::vec3(0,0,1));
    Mesh m;
    visualizer::box_to_model(&b, &m);

    std::vector<float> res;
    model_to_simple_model(&m, res);
    return res;
  }

  inline void add_model(std::vector<dfloat> &dst, const std::vector<dfloat> &src)
  {
    dst.reserve(dst.size()+src.size());
    for (int i=0;i<src.size();i++)
      dst.push_back(src[i]);
  }

  inline void add_model(std::vector<dfloat> &dst, const std::vector<float> &src)
  {
    dst.reserve(dst.size()+src.size());
    for (int i=0;i<src.size();i++)
      dst.push_back(src[i]);
  }

  inline void shift(std::vector<dfloat> &vert, dvec3 shift)
  {
    for (int i=0;i<vert.size()/FLOAT_PER_VERTEX;i++)
    {
      vert[FLOAT_PER_VERTEX*i] += shift[0];
      vert[FLOAT_PER_VERTEX*i+1] += shift[1];
      vert[FLOAT_PER_VERTEX*i+2] += shift[2];
    }
  }

  void transform(std::vector<dfloat> &vert, dmat43 mat, int floats_per_vertex = FLOAT_PER_VERTEX, int pos_start = 0, int norm_start = 3)
  {
    dmat43 norm_mat = transposedInverse3x3(mat);

    if (norm_start >= 0)
    {
      for (int i=0;i<vert.size()/floats_per_vertex;i++)
      {
        mulp(mat, vert[pos_start+floats_per_vertex*i], vert[pos_start+floats_per_vertex*i+1], vert[pos_start+floats_per_vertex*i+2]);
        mulv(norm_mat, vert[norm_start+floats_per_vertex*i], vert[norm_start+floats_per_vertex*i+1], vert[norm_start+floats_per_vertex*i+2]);

        dfloat a = vert[norm_start+floats_per_vertex*i];
        dfloat b = vert[norm_start+floats_per_vertex*i+1];
        dfloat c = vert[norm_start+floats_per_vertex*i+2];
        dfloat len = CppAD::sqrt(a*a + b*b + c*c) + 1e-18;
        vert[norm_start+floats_per_vertex*i] = a/len;
        vert[norm_start+floats_per_vertex*i+1] = b/len;
        vert[norm_start+floats_per_vertex*i+2] = c/len;
      }
    }
    else
    {
      for (int i=0;i<vert.size()/floats_per_vertex;i++)
      {
        mulp(mat, vert[pos_start+floats_per_vertex*i], vert[pos_start+floats_per_vertex*i+1], vert[pos_start+floats_per_vertex*i+2]);
      }
    }
  }

  void test_model(std::vector<dfloat> &vert, std::vector<dfloat> &params)
  {
    dvec3 shift_v{params[0], params[1], params[2]};
    dvec3 scale_v{params[3], params[4], params[5]};
    std::vector<dfloat> tri_model;
    add_model(tri_model, get_cube_expl());

    dvec3 axis{0,1,0};
    dmat43 mat = ident();
    mat = rotate(mat, axis, PI/4);
    mat = translate(mat, shift_v);
    mat = scale(mat, scale_v);
    transform(tri_model, mat);
    for (int i=0;i<4;i++)
    {
      for (int j=0;j<3;j++)
        debug("%.3f ", mat[3*i+j]);
      debugnl();
    }
    add_model(vert, tri_model);
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
      debug("%s_%.4d ", names[i % FLOAT_PER_VERTEX].c_str(), i/FLOAT_PER_VERTEX);
      for (int j = 0; j < x_n; j++)
        debug("%.2f ", jac[i * x_n + j]);
      debugnl();
    }
  }

  void dgen_test(std::vector<float> &model)
  {
    size_t x_n = 6;
    std::vector<dfloat> X(x_n);
    std::vector<int> inds;
    std::vector<dfloat> Y;

    // declare independent variables and start recording operation sequence
    logerr("gen");
    CppAD::Independent(X);

    //Cube(X[0], X[1], X[2], X[3], X[4], X[5], Y);
    test_model(Y, X);
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording
    logerr("gen_finish");
    // compute derivative using operation sequence stored in f
    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    X0[0] = 0;
    X0[1] = 0;
    X0[2] = 0;
    X0[3] = 1;
    X0[4] = 1;
    X0[5] = 1;
    jac = f.Jacobian(X0); // Jacobian for operation sequence
    res = f.Forward(0, X0);

    print_model(res);
    print_jackobian(jac, x_n, y_n);

    model = res;
  }

  bool create_model_from_block(Block &bl, ComplexModel &mod)
  {
    Model *m = new Model();
    std::vector<float> res;
    dgen::dgen_test(res);
    visualizer::simple_mesh_to_model_332(res, m);

    mod.models.push_back(m);
    mod.materials.push_back(Material(engine::textureManager->get("wood")));
    mod.update();
  }
}