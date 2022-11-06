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

  dfloat smoothmax(dfloat a, dfloat b, float alpha = 16)
  {
    return (a*exp(a*alpha) + b*exp(b*alpha))/(exp(a*alpha) + exp(b*alpha));
  }
  dfloat smoothmin(dfloat a, dfloat b, float alpha = 16)
  {
    return (a*exp(-a*alpha) + b*exp(-b*alpha))/(exp(-a*alpha) + exp(-b*alpha));
  }
  dfloat smoothclamp(dfloat a, dfloat from, dfloat to, float alpha = 16)
  {
    return smoothmin(smoothmax(a, from, alpha), to, alpha);
  }
  dfloat triangle_func(dfloat y, int n, int i)
  {
    return smoothmax(1 - abs(n * y - i), 0);
  }

  dfloat x_for_spline_y(const std::vector<dvec3> &in_spline, dfloat y, int axis_x)
  {
    dfloat sum = 0;
    for (int i = 0; i < in_spline.size(); ++i)
    {
      sum += in_spline[i][axis_x] * triangle_func(y, in_spline.size() - 1, i);
    }
    return sum;
  }

  dfloat dist(dfloat x1, dfloat y1, dfloat x2, dfloat y2)
  {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  }

  dfloat rad_by_points(const std::vector<dvec3> &in_spline, dfloat y1, dfloat y2)
  {
    return dist(x_for_spline_y(in_spline, y1, 0), y1, x_for_spline_y(in_spline, y2, 0), y2) / 2.0;
  }

  dvec3 shift_by_points(const std::vector<dvec3> &in_spline, dfloat y1, dfloat y2, dfloat thick, int x, int y)
  {
    dvec3 shift{0, 0, 0};
    shift[y] = (y1 + y2) / 2;
    shift[x] = -(x_for_spline_y(in_spline, y1, 0) + x_for_spline_y(in_spline, y2, 0) - thick) / 2;
    return shift;
  }

  dfloat sin_by_points(const std::vector<dvec3> &in_spline, dfloat y1, dfloat y2, dfloat thick)
  {
    return (x_for_spline_y(in_spline, y1, 0) - x_for_spline_y(in_spline, y2, 0) + thick / 2.0) / (2.0 * rad_by_points(in_spline, y1, y2));
  }

  std::vector<dvec3> create_spline(const std::vector<dfloat> &params, int idx, int axis_x, int axis_y, bool from_zero)
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
    for (int i=0;i<idx;i++)
    {
      dvec3 vec{0,0,0};
      vec[axis_x] = ((float)i)/(idx-1);
      vec[axis_y] = params[i];

      spline.push_back(vec);
    }

    return spline;
  }

  std::vector<dvec3> create_spline_for_handle(const std::vector<dfloat> &params, int idx, int axis_x, int axis_y)
  {
    std::vector<dvec3> spline;
    spline.reserve(1);
    dvec3 vec{0, 0, 0};
    vec[axis_y] = params[idx];
    spline.push_back(vec);
    return spline;
  }

  std::vector<dvec3> spline_rotation(const std::vector<dvec3> &in_spline, dvec3 axis, int rotations)
  {
    std::vector<dvec3> spline;
    spline.reserve(in_spline.size() * rotations);
    dmat43 rot_mat = ident();
    dfloat angle = (2 * PI) / rotations;
    rot_mat = rotate(rot_mat, axis, angle);
    for (int i = 0; i < in_spline.size(); ++i)
    {
      dvec3 vec = in_spline[i];
      for (int j = 0; j < rotations; ++j)
      {
        spline.push_back(vec);
        vec = mulp(rot_mat, vec);
      }
      spline.push_back(vec);
    }
    return spline;
  }

  std::vector<dvec3> spline_shifting(const std::vector<dvec3> &in_spline, dvec3 shift)
  {
    std::vector<dvec3> spline;
    spline.reserve(in_spline.size());
    for (int i = 0; i < in_spline.size(); ++i)
    {
      
      dvec3 vec = add(in_spline[i], shift);
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
    dvec3 d_prev = dvec3{0,0,0};
    d_prev[axis_x] = l0[axis_y];
    d_prev[axis_y] = -l0[axis_x];
    for (int i=1;i<in_spline.size()-1;i++)
    {
      dvec3 l = sub(in_spline[i+1], in_spline[i]);
      dvec3 d = dvec3{0,0,0};
      d[axis_x] = l[axis_y];
      d[axis_y] = -l[axis_x];
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

  void spline_to_model_rotate(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, int rotations)
  {
    dmat43 rot_mat = ident();
    dfloat angle = (2*PI)/rotations;
    rot_mat = rotate(rot_mat, axis, angle);
    int sp_sz = spline.size();
    int prev_size = model.size();
    model.reserve(prev_size + FLOAT_PER_VERTEX*3*2*(sp_sz-1));
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
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1), verts[i], n, dvec2{((float)(sector))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+1, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+2, verts[i-1], n, dvec2{((float)(sector))/rotations, prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+3, prev_verts[i-1], n, dvec2{((float)(sector - 1))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+4, verts[i-1], n, dvec2{((float)(sector))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+5, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        prev_len = new_len;
      }
      prev_verts = verts;
    }
  }

  void spline_to_model_part_rotate_plus_shift(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, dfloat beg_angle, dfloat part, int rotations, dvec3 shift)
  {
    dmat43 rot_mat = ident();
    dmat43 first_rot_mat = ident();
    dfloat angle = 2*part*PI/(rotations);
    rot_mat = rotate(rot_mat, axis, angle);
    dfloat ba = beg_angle + 1e-6;
    first_rot_mat = rotate(first_rot_mat, axis, ba);
    int sp_sz = spline.size();
    int prev_size = model.size();
    model.reserve(prev_size + FLOAT_PER_VERTEX*3*2*(sp_sz-1));
    std::vector<dvec3> verts = spline;
    std::vector<dvec3> prev_verts = spline;

    dfloat full_len = 1e-9;
    for (int i=1;i<sp_sz;i++)
    {
      full_len += len(sub(verts[i],verts[i-1]));
      verts[i - 1] = mulp(first_rot_mat, verts[i - 1]);
      prev_verts[i - 1] = mulp(first_rot_mat, prev_verts[i - 1]);
    }

    verts[sp_sz - 1] = mulp(first_rot_mat, verts[sp_sz - 1]);
    prev_verts[sp_sz - 1] = mulp(first_rot_mat, prev_verts[sp_sz - 1]);

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
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1), add(verts[i], shift), n, dvec2{((float)(sector))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+1, add(prev_verts[i], shift), n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+2, add(verts[i-1], shift), n, dvec2{((float)(sector))/rotations, prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+3, add(prev_verts[i-1], shift), n, dvec2{((float)(sector - 1))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+4, add(verts[i-1], shift), n, dvec2{((float)(sector))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+5, add(prev_verts[i], shift), n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        prev_len = new_len;
      }
      prev_verts = verts;
    }
  }

  void spline_to_model_part_rotate_plus_shift_with_diff_rad(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, dfloat beg_angle, dfloat angle, int rotations, dvec3 shift, int idx, const std::vector<dfloat> &params)
  {
    dmat43 rot_mat = ident();
    dfloat ba = beg_angle + 1e-6;
    angle += 1e-6;
    rot_mat = rotate(rot_mat, axis, ba);
    int sp_sz = spline.size();
    int prev_size = model.size();
    model.reserve(prev_size + FLOAT_PER_VERTEX*3*2*(sp_sz-1));
    std::vector<dvec3> verts = spline;
    std::vector<dvec3> prev_verts = spline;

    dfloat full_len = 1e-9;
    for (int i=1;i<sp_sz;i++)
    {
      full_len += len(sub(verts[i],verts[i-1]));
      prev_verts[i - 1] = mulp(rot_mat, prev_verts[i - 1]);
    }

    prev_verts[sp_sz - 1] = mulp(rot_mat, prev_verts[sp_sz - 1]);
    
    for (int sector = 1; sector <= rotations; sector++)
    {
      rot_mat = rotate(rot_mat, axis, angle);
      dvec3 vert;
      dfloat coef = params[idx + sector] / params[idx];//now radius / first radius
      for (int i=0;i<sp_sz;i++)
      {
        vert = {spline[i][0], spline[i][1] * coef, spline[i][2]};
        verts[i] = mulp(rot_mat, vert);
      }
      dfloat prev_len = 0;
      dfloat new_len = 0;
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + len(sub(verts[i],verts[i-1]));
        dvec3 v1 = sub(prev_verts[i], verts[i]);
        dvec3 v2 = sub(verts[i-1], verts[i]); 
        dvec3 n = normalize(cross(v1, v2));
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1), add(verts[i], shift), n, dvec2{((float)(sector))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+1, add(prev_verts[i], shift), n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+2, add(verts[i-1], shift), n, dvec2{((float)(sector))/rotations, prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+3, add(prev_verts[i-1], shift), n, dvec2{((float)(sector - 1))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+4, add(verts[i-1], shift), n, dvec2{((float)(sector))/rotations,  prev_len/full_len}); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+5, add(prev_verts[i], shift), n, dvec2{((float)(sector - 1))/rotations, new_len/full_len}); 
        prev_len = new_len;
      }
      prev_verts = verts;
    }
  }

  void create_cup(const std::vector<dfloat> &params, std::vector<dfloat> &vert)
  {
    std::vector<dvec3> spline = create_spline(params, 9, 1, 0, true);
    dmat43 sc = scale(ident(), dvec3{0.09,0.9,0.09});
    transform(spline, sc);
    //spline = spline_make_smoother(spline, 4, 1, -1, 1, 0);
    if (params[10] > 0.5)
    {
      int handle_param_idx = 11;
      std::vector<dvec3> spline1 = create_spline_for_handle(params, handle_param_idx, 0, 1);
      dfloat thick = smoothmin(params[handle_param_idx], 0.02, 8);
      dfloat start_pos = params[handle_param_idx+1] - params[handle_param_idx];
      dfloat end_pos = params[handle_param_idx+1] + params[handle_param_idx+2] + params[handle_param_idx] + thick;
      spline1 = spline_rotation(spline1, dvec3{1, 0, 0}, 8);
      spline1 = spline_shifting(spline1, dvec3{0, rad_by_points(spline, start_pos, end_pos), 0});
      dfloat sin_p = smoothmin(sin_by_points(spline, end_pos, (start_pos + end_pos) / 2.0, thick), 0.98, 8);
      spline_to_model_part_rotate_plus_shift(vert, spline1, dvec3{0, 0, 1}, asin(sin_p), 0.5, 8, shift_by_points(spline, start_pos, end_pos, thick, 0, 1));
    }
    spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, spline, dvec3{0,1,0},16);
    dmat43 sc2 = scale(ident(), dvec3{1,params[9],1});
    transform(vert, sc2);
  }

  void create_cup_2(const std::vector<dfloat> &params, std::vector<dfloat> &vert)
  {
    std::vector<dvec3> spline = create_spline(params, 9, 1, 0, true);
    dmat43 sc = scale(ident(), dvec3{0.09,0.9,0.09});
    transform(spline, sc);
    if (params[10] > 0.5)
    {
      int handle_param_idx = 11;
      std::vector<dvec3> spline1 = create_spline_for_handle(params, handle_param_idx, 0, 1);//thick
      spline1 = spline_rotation(spline1, dvec3{1, 0, 0}, 8);
      spline1 = spline_shifting(spline1, dvec3{0, params[handle_param_idx + 5], 0});//first radius
      dvec3 center = {params[handle_param_idx + 1], params[handle_param_idx + 2], 0};//center coords
      dfloat alpha = params[handle_param_idx + 3];//step angle
      dfloat start = params[handle_param_idx + 4];//beg angle
      spline_to_model_part_rotate_plus_shift_with_diff_rad(vert, spline1, dvec3{0, 0, 1}, start, alpha, 8, center, handle_param_idx + 5, params);
    }
    spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, spline, dvec3{0,1,0},16);
    dmat43 sc2 = scale(ident(), dvec3{1,params[9],1});
    transform(vert, sc2);
  }

  dfloat parameters_limits_reg(const std::vector<dfloat> &params, const std::vector<float> &params_min, const std::vector<float> &params_max,
                               float edge_size)
  {
    dfloat res = 0;
    for (int i = 0; i < params.size(); i++)
    {
      res += smoothmax((params_min[i] + edge_size - params[i])/edge_size, 0.0f, 8);
      res += smoothmax((params[i] - (params_max[i] - edge_size))/edge_size, 0.0f, 8);
    }
    res = smoothmax(smoothmax(smoothmax(res, 0, 8), 0, 8), 0, 8);
    return res;
  }
  dfloat parameters_cup_reg(const std::vector<dfloat> &params)
  {
    int spline_offsets_cnt = 9;
    dfloat res = 0;
    for (int i = 1; i < spline_offsets_cnt; i++)
    {
      res += smoothmax(params[i-1] - params[i] - 0.01f, 0, 8);
    }
    res = smoothmax(smoothmax(smoothmax(res, 0, 8), 0, 8), 0, 8);
    return res;
  }

  void transform_by_scene_parameters(std::vector<dgen::dfloat> &params, int offset, std::vector<dgen::dfloat> &model)
  {
    dmat43 rot = rotate(ident(), dvec3{1,0,0}, params[offset]);
    rot = rotate(rot, dvec3{0,1,0}, params[offset+1]);
    rot = rotate(rot, dvec3{0,0,1}, params[offset+2]);
    rot = translate(rot, dvec3{params[offset+3], params[offset+4], params[offset+5]});
    transform(model, rot);
  }

  bool check_stability(generator_func func, const std::vector<float> &params, int iterations)
  {
    float eps = 1e-4;
    std::vector<float> model_ref;
    std::vector<float> jac_ref;
    dgen_test_internal(model_ref, func, params, params, &jac_ref);
    int x_n = params.size();
    int y_n = model_ref.size();
    debug("Checking stability of differential procedural model\n");
    debug("Model has %d params and %d vertices\n", x_n, y_n/FLOAT_PER_VERTEX);

    std::vector<int> failed_params;
    for (int param_n = 0; param_n < x_n; param_n++)
    {
      int failed_tests = 0;
      for (int i=0;i<iterations;i++)
      {
        std::vector<float> par = params;
        par[param_n] *= urand(0.5, 2);
        std::vector<float> model;
        std::vector<float> jac;
        bool model_created = false;
        try
        {
          dgen_test_internal(model, func, par, params, &jac);
          model_created = true;
        }
        catch(const std::exception& e)
        {
          std::cerr << e.what() << '\n';
        }
        if (!model_created)
        {
          debug("Test %d failed. Generator crashed", i);
        }
        else if (model.size() != model_ref.size())
        {
          debug("Test %d failed. Model has wrong number of vertices (%d)\n", i, model.size()/FLOAT_PER_VERTEX);
          failed_tests++;
        }
        else
        {
          int diff_mod = 0;
          for (int j=0;j<model.size();j++)
          {
            if (abs(model[j] - model_ref[j])/(abs(model[j] + model_ref[j]) + eps) > eps)
              diff_mod++;
          }
          int diff_jac = 0;
          for (int j=0;j<jac.size();j++)
          {
            if (abs(jac[j] - jac_ref[j])/(abs(jac[j] + jac_ref[j]) + eps) > eps)
              diff_jac++;
          }
          if (diff_jac > 0 || diff_mod > 0)
            failed_tests++;
          if (diff_mod > 0)
          {
            debug("Test %d failed. Model has difference in %d values with the base one\n", i, diff_mod);
          }
          if (diff_jac > 0)
          {
            debug("Test %d failed. Jacobian has difference in %d values with the base one\n", i, diff_jac);
          }
        }
      }

      if (failed_tests == 0)
      {
        debug("Param %d. Stability check PASSED\n", param_n);
      }
      else
      {
        debug("Param %d. Stability check FAILED (%d/%d fails)\n", param_n, failed_tests, iterations);
        failed_params.push_back(param_n);
      }
    }

    if (failed_params.empty())
    {
      debug("Stability check PASSED\n");
      return true;
    }
    else
    {
      debug("Stability check FAILED\n");
      debug("%d unstable parameters: {", failed_params.size());
      for (int i=0;i<failed_params.size();i++)
      {
        debug("%d", failed_params[i]);
        if (i != failed_params.size() - 1)
          debug(", ");
      }
      debug("}\n");

      return false;
    }
  }

  bool check_robustness(generator_func func, const std::vector<float> &params_min, const std::vector<float> &params_max, int iterations)
  {
    for (int i=0;i<iterations;i++)
    {
      std::vector<float> X0;
      std::vector<float> model;
      for (int j=0;j<params_min.size();j++)
      {
        float rnd = urand(0,1);
        X0.push_back(rnd*params_min[j] + (1-rnd)*params_max[j]);
      }
      dgen_test_internal(model, func, X0, X0);
    }
    return true;
  }

  void dgen_test(std::vector<float> &params, std::vector<float> &model)
  {
    dgen_test_internal(model, create_cup, params, params);
  }
  void dgen_test_internal(std::vector<float> &model, generator_func func, const std::vector<float> &check_params, 
                          const std::vector<float> &params, std::vector<float> *jacobian )
  {
    assert(check_params.size() > 0);
    assert(check_params.size() == params.size());
    assert(model.empty());

    size_t x_n = check_params.size();
    std::vector<dfloat> X(x_n);
    std::vector<dfloat> Y;
    for (int i=0;i<x_n;i++)
      X[i] = check_params[i];

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X);
    func(X, Y);
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y); // store operation sequence in f: X -> Y and stop recording

    // compute derivative using operation sequence stored in f
    std::vector<float> jac(y_n * x_n); // Jacobian of f (m by n matrix)
    std::vector<float> res(y_n); 
    std::vector<float> X0(x_n);        // domain space vector
    for (int i=0;i<x_n;i++)
      X0[i] = params[i];
  
    jac = f.Jacobian(X0); // Jacobian for operation sequence
    res = f.Forward(0, X0);

    //print_model(res);
    //print_jackobian(jac, x_n, y_n);

    model = res;
    if (jacobian)
      *jacobian = jac;
  }

  bool create_model_from_block(Block &bl, ComplexModel &mod)
  {
    Model *m = new Model();
        std::vector<float> X0{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                          1, 1, 
                          0.08, 0.25, 0.5};//handle params
    std::vector<float> res;
    dgen::dgen_test(X0, res);
    visualizer::simple_mesh_to_model_332(res, m);

    mod.models.push_back(m);
    mod.materials.push_back(Material(engine::textureManager->get("porcelain")));
    mod.update();
    return true;
  }
}