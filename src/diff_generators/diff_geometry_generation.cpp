#include "diff_geometry_generation.h"
#include "common_utils/utility.h"
#include "graphics_utils/modeling.h"
#include <cppad/cppad.hpp>
#include "tinyEngine/engine.h"
#include "common_utils/interpolation.h"
namespace dgen
{
  inline void add_vertex(std::vector<dfloat> &vert, int n, const dvec3 &pos, const dvec3 &norm, const dvec2 &tc)
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

  void transform(std::vector<dfloat> &vert, dmat43 mat, int floats_per_vertex, int pos_start, int norm_start)
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

  void transform(std::vector<dvec3> &verts, dmat43 mat)
  {
    for (auto &vert : verts)
    {
      vert = mulp(mat, vert); 
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

  void print_jackobian(const std::vector<float> &jac, int x_n, int y_n, int lines)
  {
    std::string names[FLOAT_PER_VERTEX] = {" pos_x", " pos_y", " pos_z", 
                                           "norm_x", "norm_y", "norm_z",
                                           "  tc_x", "  tc_y"};
    debug("Jacobian    ");
    for (int j = 0; j < x_n; j++)
      debug("x_%.2d ", j);
    debugnl();
    for (int i = 0; i < MIN(y_n, lines); i++)
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
  std::vector<dvec3> create_spline(const std::vector<dfloat> &params, int axis_x, int axis_y, bool from_zero)
  {
    //if you have n params and axis_x = 0, axis_y = 2 then it will create n points (vec3) with formula point[i] = (i/(n-1), 0, params[i])
    std::vector<dvec3> spline;
    spline.reserve((params.size() + from_zero));
    if (from_zero)
    {
      dvec3 vec{0,0,0};
      vec[axis_x] = -1e-4;
      spline.push_back(vec);
    }
    for (int i=0;i<params.size();i++)
    {
      dvec3 vec{0,0,0};
      vec[axis_x] = ((float)i)/(params.size()-1);
      vec[axis_y] = params[i];

      spline.push_back(vec);
    }

    return spline;
  }

  std::vector<dvec3> spline_to_closed_curve_thickness(const std::vector<dvec3> &in_spline, dfloat thickness, int axis_x, int axis_y)
  {
    std::vector<dvec3> spline = in_spline;
    std::vector<dvec3> offset_dirs;
    dvec3 x = dvec3{0,0,0};
    x[axis_x] = 1;
    dvec3 y = dvec3{0,0,0};
    y[axis_y] = -1;
    offset_dirs.push_back(add(in_spline[0],mul(thickness,x)));

    dvec3 l0 = sub(in_spline[1], in_spline[0]);
    dvec3 n0 = cross(l0, x);
    dvec3 d_prev = cross(n0, l0);

    for (int i=1;i<in_spline.size()-1;i++)
    {
      dvec3 l = sub(in_spline[i+1], in_spline[i]);
      dvec3 n = cross(l, x);
      dvec3 d = cross(n, l);
      offset_dirs.push_back(add(in_spline[i],mul(thickness,normalize(add(d_prev, d))))); //pos + thickness*dir
      d_prev = d;
    }

    offset_dirs.push_back(add(in_spline[in_spline.size()-1],mul(thickness,y)));

    for (int i=in_spline.size()-1; i>=0; i--)
    {
      spline.push_back(offset_dirs[i]);
    }

    return spline;
  }

  std::vector<dvec3> spline_to_closed_curve_mirror(const std::vector<dvec3> &in_spline, int axis_x, int axis_y)
  {
    std::vector<dvec3> spline = in_spline;
    for (int i=in_spline.size()-1; i>=0; i--)
    {
      dvec3 vec_m = in_spline[i];
      vec_m[axis_y] *= -1;
      spline.push_back(vec_m);
    }
    return spline;
  }

  std::vector<dvec3> spline_make_smoother(const std::vector<dvec3> &in_spline, int detail_q, int index_from, int index_to, int axis_x, int axis_y)
  {
    index_from = index_from < 0 ? 0 : index_from;
    index_to = index_to < 0 ? in_spline.size() : index_to;
    if (detail_q <= 1)
      return in_spline;
    std::vector<dfloat> xs;
    std::vector<dfloat> ys;
    for (int i = index_from; i < index_to; i++)
    {
      xs.push_back(in_spline[i][axis_x]);
      ys.push_back(in_spline[i][axis_y]);
    }
    std::vector<interpolation::Spline<dfloat>> splines = interpolation::spline<dfloat>(xs, ys);
    std::vector<dvec3> sp;
    dfloat step = 1.0f/(detail_q+1);

    for (int i = 0; i < index_from; i++)
      sp.push_back(in_spline[i]);
    for (int s_n = 0; s_n < splines.size(); s_n++)
    {
      for (int i=0;i<detail_q;i++)
      {
        dfloat pos = splines[s_n].x + i*step*(xs[s_n+1] - xs[s_n]);
        dfloat val = splines[s_n].get(pos);
        dvec3 vec{0,0,0};
        vec[axis_x] = pos;
        vec[axis_y] = val;

        sp.push_back(vec);
      }
    }
    for (int i = index_to; i < in_spline.size(); i++)
      sp.push_back(in_spline[i]);

    return sp;
  }

  void spline_to_model_rotate(std::vector<dfloat> &model, int &model_size, const std::vector<dvec3> &spline, dvec3 axis, int rotations)
  {
    dmat43 rot_mat = ident();
    dfloat angle = (2*PI)/rotations;
    rot_mat = rotate(rot_mat, axis, angle);
    int sp_sz = spline.size();
    int offset = model_size / FLOAT_PER_VERTEX;
    int part_size = FLOAT_PER_VERTEX*3*2*(sp_sz-1);
    model_size += part_size;
    model.reserve(model_size);
    std::vector<dvec3> verts = spline;
    std::vector<dvec3> prev_verts = spline;

    dfloat full_len = 1e-9;
    for (int i=1;i<sp_sz;i++)
    {
      full_len += len(sub(verts[i],verts[i-1]));
    }

    for (int sector = 1; sector <= rotations; sector++)
    {
      for (int i=0;i<sp_sz;i++)
      {
        verts[i] = mulp(rot_mat, verts[i]);
      }
      dfloat prev_len = 0;
      dfloat new_len = 0;
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + len(sub(verts[i],verts[i-1]));
        dvec3 v1 = sub(prev_verts[i], verts[i]);
        dvec3 v2 = sub(verts[i-1], verts[i]); 
        dvec3 n = normalize(cross(v1, v2));
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1), verts[i], n, dvec2{((float)(sector))/rotations, new_len/full_len}); 
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1)+1, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1)+2, verts[i-1], n, dvec2{((float)(sector))/rotations, prev_len/full_len}); 
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1)+3, prev_verts[i-1], n, dvec2{((float)(sector - 1))/rotations,  prev_len/full_len}); 
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1)+4, verts[i-1], n, dvec2{((float)(sector))/rotations,  prev_len/full_len}); 
        add_vertex(model, offset + 6*((sp_sz-1)*(sector-1) + i - 1)+5, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        prev_len = new_len;
      }
      prev_verts = verts;
    }
  }

  void test_spline(std::vector<dfloat> &vert, std::vector<dfloat> &params)
  {
    int model_size = 0;
    std::vector<dvec3> spline = create_spline(params, 1, 0, true);
    dmat43 sc = scale(ident(), dvec3{0.1,1,0.1});
    transform(spline, sc);
    spline = spline_make_smoother(spline, 4, 1, -1, 1, 0);
    spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, model_size, spline, dvec3{0,1,0},32);
  }

  void create_cup(std::vector<dfloat> &params, std::vector<dfloat> &vert)
  {
    int model_size = 0;
    std::vector<dvec3> spline = create_spline(params, 1, 0, true);
    dmat43 sc = scale(ident(), dvec3{0.1,1,0.1});
    transform(spline, sc);
    //spline = spline_make_smoother(spline, 4, 1, -1, 1, 0);
    //spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, model_size, spline, dvec3{0,1,0},8);
  }

  void dgen_test(std::vector<float> &model)
  {
    size_t x_n = 9;
    std::vector<dfloat> X(x_n);
    std::vector<int> inds;
    std::vector<dfloat> Y;

    // declare independent variables and start recording operation sequence
    logerr("gen");
    CppAD::Independent(X);

    //Cube(X[0], X[1], X[2], X[3], X[4], X[5], Y);
    test_spline(Y, X);
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording
    logerr("gen_finish");
    // compute derivative using operation sequence stored in f
    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    
    X0[0] = 4 - 1.45;
    X0[1] = 4 - 1.0;
    X0[2] = 4 - 0.65;
    X0[3] = 4 - 0.45;
    X0[4] = 4 - 0.25;
    X0[5] = 4 - 0.18;
    X0[6] = 4 - 0.1;
    X0[7] = 4 - 0.05;
    X0[8] = 4 - 0;
    //X0[9] = 0.3 + 0.81;
    //X0[10] = 0.3 + 1.0;
    //X0[11] = 0.3 + 1.21;
    //X0[12] = 0.3 + 1.44;
    //X0[13] = 0.3 + 1.69;
  
    jac = f.Jacobian(X0); // Jacobian for operation sequence
    res = f.Forward(0, X0);

    print_model(res);
    //print_jackobian(jac, x_n, y_n);

    model = res;
  }

  bool create_model_from_block(Block &bl, ComplexModel &mod)
  {
    Model *m = new Model();
    std::vector<float> res;
    dgen::dgen_test(res);
    visualizer::simple_mesh_to_model_332(res, m);

    mod.models.push_back(m);
    mod.materials.push_back(Material(engine::textureManager->get("porcelain")));
    mod.update();
    return true;
  }
}