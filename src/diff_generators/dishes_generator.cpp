#include "dishes_generator.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  template<typename float_type>
  float_type triangle_func(float_type y, int n, int i)
  {
    if constexpr (std::is_same<float_type, dfloat>::value)
    {
      return d_max(1 - abs(n * y - i), 0);
    }
    return MAX(1 - abs(n * y - i), 0);
  }
  template dfloat triangle_func<dfloat>(dfloat y, int n, int i);

  template<typename float_type>
  float_type x_for_spline_y(const std::vector<g_vec3<float_type>> &in_spline, float_type y, int start_id)
  {
    float_type sum = 0;
    if (y < in_spline[0][1])
      return in_spline[0][0];
    if (y > in_spline.back()[1])
      return in_spline.back()[0];

    for (int i = start_id + 1; i < in_spline.size(); ++i)
    {
      if (y >= in_spline[i-1][1] && y <= in_spline[i][1])
      {
        float_type d = in_spline[i][1] - in_spline[i-1][1];
        float_type d1 = ((y - in_spline[i-1][1])/d);
        //std::cerr<<i<<")"<<in_spline[i][0]<<" "<<in_spline[i][1]<<" "<<y<<" "<<start_id<<" "<<sum<<"\n";
        return (1-d1)*in_spline[i-1][0] + d1*in_spline[i][0];
      }
      //sum += in_spline[i][0] * triangle_func(y, in_spline.size() - 1, i);
      //std::cerr<<in_spline[i][0]<<" "<<in_spline[i][1]<<" "<<y<<" "<<start_id<<" "<<sum<<"\n";
    }
    return in_spline.back()[0];
    //std::cerr<<"\n";
    //return sum;
  }
  template dfloat x_for_spline_y<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, dfloat y, int start_id);

  template<typename float_type>
  float_type dist(float_type x1, float_type y1, float_type x2, float_type y2)
  {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  }
  template dfloat dist<dfloat>(dfloat x1, dfloat y1, dfloat x2, dfloat y2);

  template<typename float_type>
  float_type rad_by_points(const std::vector<g_vec3<float_type>> &in_spline, float_type y1, float_type y2, int start_id)
  {
    return dist(x_for_spline_y(in_spline, y1, start_id), y1, x_for_spline_y(in_spline, y2, start_id), y2) / 2.0;
  }
  template dfloat rad_by_points<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, dfloat y1, dfloat y2, int start_id);

  template<typename float_type>
  g_vec3<float_type> shift_by_points(const std::vector<g_vec3<float_type>> &in_spline, float_type center_y, float_type thick, int x, int y)
  {
    g_vec3<float_type> shift{0, 0, 0};
    shift[y] = center_y;
    shift[x] = -x_for_spline_y(in_spline, center_y, 0) + thick;
    return shift;
  }
  template g_vec3<dfloat> shift_by_points<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, dfloat center_y, dfloat thick, int x, int y);

  template<typename float_type>
  float_type sin_by_points(const std::vector<g_vec3<float_type>> &in_spline, float_type y1, float_type y2, float_type thick, int start_id)
  {
    return (x_for_spline_y(in_spline, y1, start_id) - x_for_spline_y(in_spline, y2, start_id) + thick / 2.0) / (2.0 * rad_by_points(in_spline, y1, y2, start_id));
  }
  template dfloat sin_by_points<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, dfloat y1, dfloat y2, dfloat thick, int start_id);

  template <typename float_type>
  std::vector<g_vec3<float_type>> create_spline(const std::vector<float_type> &params, int idx, int axis_x, int axis_y, bool from_zero, int detail_q)
  {
    //if you have n params and axis_x = 0, axis_y = 2 then it will create n points (vec3) with formula point[i] = (i/(n-1), 0, params[i])
    std::vector<g_vec3<float_type>> spline;
    spline.reserve((params.size() + from_zero));
    if (from_zero)
    {
      float_type smooth_y;
      if constexpr (std::is_same<float_type, dfloat>::value)
      {
        smooth_y = d_min(0.5, 0.9*params[0]);
      }
      else
      {
        smooth_y = MIN(0.5, 0.9*params[0]);
      }
      float_type smooth_x = 0.05;
      //center of cup bottom
      g_vec3<float_type> vec{0,0,0};
      vec[axis_x] = -smooth_x -1e-4;
      spline.push_back(vec);

      //smoothness on the bottom/side transition
      int steps = 2 + 2*detail_q;
      for (int i=0;i<steps;i++)
      {
        g_vec3<float_type> vec{0,0,0};
        vec[axis_x] = -smooth_x*cos((PI*i)/(2*steps));
        vec[axis_y] = params[0] - smooth_y*(1 - sin((PI*i)/(2*steps)));
        spline.push_back(vec);
      }
    }
    for (int i=0;i<idx;i++)
    {
      g_vec3<float_type> vec{0,0,0};
      vec[axis_x] = ((double)i)/(idx-1);
      vec[axis_y] = params[i];

      spline.push_back(vec);
    }

    return spline;
  }
  template std::vector<g_vec3<dfloat>> create_spline<dfloat>(const std::vector<dfloat> &params, int idx, int axis_x, int axis_y, bool from_zero, int detail_q);

  template<typename float_type>
  std::vector<g_vec3<float_type>> create_spline_for_handle(int axis_x, int axis_y)
  {
    std::vector<g_vec3<float_type>> spline;
    spline.reserve(1);
    g_vec3<float_type> vec{0, 0, 0};
    vec[axis_y] = 1;
    spline.push_back(vec);
    return spline;
  }
  template std::vector<g_vec3<dfloat>> create_spline_for_handle<dfloat>(int axis_x, int axis_y);

  template<typename float_type>
  std::vector<g_vec3<float_type>> spline_rotation(const std::vector<g_vec3<float_type>> &in_spline, g_vec3<float_type> axis, int rotations)
  {
    std::vector<g_vec3<float_type>> spline;
    spline.reserve(in_spline.size() * rotations);
    g_mat43<float_type> rot_mat = ident<float_type>();
    float_type angle = (2 * PI) / rotations;
    rot_mat = rotate(rot_mat, axis, angle);
    for (int i = 0; i < in_spline.size(); ++i)
    {
      g_vec3<float_type> vec = in_spline[i];
      for (int j = 0; j < rotations; ++j)
      {
        spline.push_back(vec);
        vec = mulp(rot_mat, vec);
      }
      spline.push_back(vec);
    }
    return spline;
  }
  template std::vector<g_vec3<dfloat>> spline_rotation<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, g_vec3<dfloat> axis, int rotations);

  template<typename float_type>
  std::vector<g_vec3<float_type>> spline_to_closed_curve_thickness(const std::vector<g_vec3<float_type>> &in_spline, float_type thickness, int axis_x, int axis_y)
  {
    std::vector<g_vec3<float_type>> spline = in_spline;
    std::vector<g_vec3<float_type>> offset_dirs;
    g_vec3<float_type> x = g_vec3<float_type>{0,0,0};
    x[axis_x] = 1;
    g_vec3<float_type> y = g_vec3<float_type>{0,0,0};
    y[axis_y] = -1;
    offset_dirs.push_back(in_spline[0] + thickness*x);

    g_vec3<float_type> l0 = in_spline[1] - in_spline[0];
    g_vec3<float_type> d_prev = g_vec3<float_type>{0,0,0};
    d_prev[axis_x] = l0[axis_y];
    d_prev[axis_y] = -l0[axis_x];
    for (int i=1;i<in_spline.size()-1;i++)
    {
      g_vec3<float_type> l = in_spline[i+1] - in_spline[i];
      g_vec3<float_type> d = g_vec3<float_type>{0,0,0};
      d[axis_x] = l[axis_y];
      d[axis_y] = -l[axis_x];
      offset_dirs.push_back(in_spline[i] + thickness*normalize(d_prev + d)); //pos + thickness*dir
      d_prev = d;
    }
    offset_dirs.push_back(in_spline[in_spline.size()-1] + thickness*y);

    //smooth transition on top
    float_type d = spline.back()[axis_y] - offset_dirs.back()[axis_y];
    {
      g_vec3<float_type> v = spline.back();
      v[axis_x] += 0.25*d;
      v[axis_y] -= 0.25*d;
      spline.push_back(v);
    }
    {
      g_vec3<float_type> v = spline.back();
      v[axis_x] += 0.25*d;
      v[axis_y] -= 0.75*d;
      spline.push_back(v);
    }
    for (int i=in_spline.size()-1; i>=0; i--)
    {
      spline.push_back(offset_dirs[i]);
    }

    return spline;
  }
  template std::vector<g_vec3<dfloat>> spline_to_closed_curve_thickness<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, dfloat thickness, int axis_x, int axis_y);

  template<typename float_type>
  std::vector<g_vec3<float_type>> spline_to_closed_curve_mirror(const std::vector<g_vec3<float_type>> &in_spline, int axis_x, int axis_y)
  {
    std::vector<g_vec3<float_type>> spline = in_spline;
    for (int i=in_spline.size()-1; i>=0; i--)
    {
      g_vec3<float_type> vec_m = in_spline[i];
      vec_m[axis_y] *= -1;
      spline.push_back(vec_m);
    }
    return spline;
  }
  template std::vector<g_vec3<dfloat>> spline_to_closed_curve_mirror<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, int axis_x, int axis_y);

  template<typename float_type>
  float_type lagrange_3_dots(float_type x, float_type x1, float_type y1, float_type x2, float_type y2, float_type x3, float_type y3)
  {
    float_type y =  y1*((x - x2)/(x1 - x2))*((x - x3)/(x1 - x3)) +
           y2*((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3)) +
           y3*((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2));
    //std::cerr<<x<<"--"<<y<<"--"<<x1<<"--"<<x2<<"--"<<x3<<"--"<<((x - x2)/(x1 - x2))*((x - x3)/(x1 - x2))<<"--"<<((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3))<<"--"<<((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2))<<"\n";
    return y;
  }
  template dfloat lagrange_3_dots<dfloat>(dfloat x, dfloat x1, dfloat y1, dfloat x2, dfloat y2, dfloat x3, dfloat y3);

  template<typename float_type>
  std::vector<g_vec3<float_type>> spline_make_smoother(const std::vector<g_vec3<float_type>> &in_spline, int detail_q, int index_from, int index_to, int axis_x, int axis_y)
  {
    index_from = index_from < 0 ? 0 : index_from;
    index_to = index_to < 0 ? in_spline.size() : index_to;
    if (detail_q <= 1)
      return in_spline;
    if (detail_q == 2 && (in_spline.size() == index_to))
    {
      std::vector<g_vec3<float_type>> sp = std::vector<g_vec3<float_type>>(index_from + 2*(index_to - index_from) - 1, g_vec3<float_type>{0,0,0});
      for (int i = 0; i < index_from; i++)
        sp[i] = in_spline[i];
      for (int i = index_from; i < index_to; i++)
        sp[index_from + 2*(i - index_from)] = in_spline[i];
      for (int i = index_from+1; i < index_to-1; i++)
      {
        float_type new_x1 = (in_spline[i-1][axis_x] + in_spline[i][axis_x])/2;
        float_type new_x2 = (in_spline[i][axis_x] + in_spline[i+1][axis_x])/2;
        sp[index_from + 2*(i - index_from) - 1][axis_x] = new_x1;
        sp[index_from + 2*(i - index_from) - 1][axis_y] +=0.5*lagrange_3_dots(new_x1, 
                                                                              in_spline[i-1][axis_x], in_spline[i-1][axis_y],
                                                                              in_spline[i  ][axis_x], in_spline[i  ][axis_y],
                                                                              in_spline[i+1][axis_x], in_spline[i+1][axis_y]);
        sp[index_from + 2*(i - index_from) + 1][axis_x] = new_x2;
        sp[index_from + 2*(i - index_from) + 1][axis_y] +=0.5*lagrange_3_dots(new_x2, 
                                                                              in_spline[i-1][axis_x], in_spline[i-1][axis_y],
                                                                              in_spline[i  ][axis_x], in_spline[i  ][axis_y],
                                                                              in_spline[i+1][axis_x], in_spline[i+1][axis_y]);
      }
      sp[index_from + 1][axis_y] *= 2;
      sp[index_from + 2*(index_to-2 - index_from) + 1][axis_y] *= 2;
      //for (int i = 0; i < sp.size(); i++)
      //  std::cerr<<"["<<i<<"]"<<sp[i][0]<<" "<<sp[i][1]<<" "<<sp[i][2]<<"\n";
      return sp;
    }
    std::vector<float_type> xs;
    std::vector<float_type> ys;
    for (int i = index_from; i < index_to; i++)
    {
      xs.push_back(in_spline[i][axis_x]);
      ys.push_back(in_spline[i][axis_y]);
    }
    std::vector<interpolation::Spline<float_type>> splines = interpolation::spline<float_type>(xs, ys);
    std::vector<g_vec3<float_type>> sp;
    float_type step = 1.0f/(detail_q+1);

    for (int i = 0; i < index_from; i++)
      sp.push_back(in_spline[i]);
    for (int s_n = 0; s_n < splines.size(); s_n++)
    {
      for (int i=0;i<detail_q;i++)
      {
        float_type pos = splines[s_n].x + i*step*(xs[s_n+1] - xs[s_n]);
        float_type val = splines[s_n].get(pos);
        g_vec3<float_type> vec{0,0,0};
        vec[axis_x] = pos;
        vec[axis_y] = val;

        sp.push_back(vec);
      }
    }
    for (int i = index_to; i < in_spline.size(); i++)
      sp.push_back(in_spline[i]);

    return sp;
  }
  template std::vector<g_vec3<dfloat>> spline_make_smoother<dfloat>(const std::vector<g_vec3<dfloat>> &in_spline, int detail_q, int index_from, int index_to, int axis_x, int axis_y);

  template<typename float_type>
  void spline_to_model_rotate(std::vector<float_type> &model, const std::vector<g_vec3<float_type>> &spline, g_vec3<float_type> axis, int rotations, bool only_pos)
  {
    g_mat43<float_type> rot_mat = ident<float_type>();
    float_type angle = (2*PI)/rotations;
    rot_mat = rotate(rot_mat, axis, angle);
    int sp_sz = spline.size();
    int prev_size = model.size()/FLOAT_PER_VERTEX;
    model.reserve(prev_size + FLOAT_PER_VERTEX*3*2*(sp_sz-1));
    std::vector<g_vec3<float_type>> verts = spline;
    std::vector<g_vec3<float_type>> prev_verts = spline;

    float_type full_len = 1e-9;
    for (int i=1;i<sp_sz;i++)
    {
      full_len += length(verts[i]-verts[i-1]);
    }

    std::vector<g_vec3<float_type>> positions;//size 4*(sp_sz-1)*rotations
    positions.reserve(4*(sp_sz-1)*rotations);
    std::vector<g_vec3<float_type>> normals;//size (sp_sz-1)*rotations
    normals.reserve((sp_sz-1)*rotations);
    std::vector<g_vec2<float_type>> tc;//size 4*(sp_sz-1)*rotations
    tc.reserve(4*(sp_sz-1)*rotations);
    std::vector<int> indices;//size (sp_sz-1)*rotations
    indices.reserve((sp_sz-1)*rotations);

    for (int sector = 1; sector <= rotations; sector++)
    {
      for (int i=0;i<sp_sz;i++)
      {
        verts[i] = mulp(rot_mat, verts[i]);
      }
      float_type prev_len = 0;
      float_type new_len = 0;
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + length(verts[i] - verts[i-1]);
        g_vec3<float_type> v1 = prev_verts[i-1] - verts[i];
        g_vec3<float_type> v2 = verts[i-1] - prev_verts[i]; 
        g_vec3<float_type> n = normalize(cross(v1, v2));

        positions.push_back(verts[i]);
        positions.push_back(prev_verts[i]);
        positions.push_back(verts[i-1]);
        positions.push_back(prev_verts[i-1]);

        normals.push_back(n);

        tc.push_back(g_vec2<float_type>{((float)(sector))/rotations, 0.75*new_len/full_len});
        tc.push_back(g_vec2<float_type>{((float)(sector - 1))/rotations, 0.75*new_len/full_len});
        tc.push_back(g_vec2<float_type>{((float)(sector))/rotations, 0.75*prev_len/full_len});
        tc.push_back(g_vec2<float_type>{((float)(sector - 1))/rotations,  0.75*prev_len/full_len});

        indices.push_back(prev_size + 6*((sp_sz-1)*(sector-1) + i - 1));
        
        prev_len = new_len;
      }
      prev_verts = verts;
    }

    for (int i=0;i<(sp_sz-1)*rotations;i++)
    {
      int sp_n = i % (sp_sz-1) + 1;
      int rot_n = i / (sp_sz-1) + 1;

      g_vec3<float_type> n = normals[i];
      g_vec3<float_type> n01 = (sp_n == sp_sz-1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n_11 = (sp_n == sp_sz-1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n11 = (sp_n == sp_sz-1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n_10 = normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 1];
      g_vec3<float_type> n10 = normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 1];
      g_vec3<float_type> n0_1 = (sp_n == 1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n - 2];
      g_vec3<float_type> n_1_1 = (sp_n == 1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 2];
      g_vec3<float_type> n1_1 = (sp_n == 1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 2];

      g_vec3<float_type> n0 = normalize(n + n01 + n11 + n10);
      g_vec3<float_type> n1 = normalize(n + n01 + n_11 + n_10);
      g_vec3<float_type> n2 = normalize(n + n10 + n1_1 + n0_1);
      g_vec3<float_type> n3 = normalize(n + n0_1 + n_1_1 + n_10);
      add_vertex<float_type>(model, indices[i]+0, positions[4*i+0], n0, tc[4*i+0], only_pos); 
      add_vertex<float_type>(model, indices[i]+1, positions[4*i+1], n1, tc[4*i+1], only_pos); 
      add_vertex<float_type>(model, indices[i]+2, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex<float_type>(model, indices[i]+3, positions[4*i+3], n3, tc[4*i+3], only_pos); 
      add_vertex<float_type>(model, indices[i]+4, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex<float_type>(model, indices[i]+5, positions[4*i+1], n1, tc[4*i+1], only_pos); 
    }
  }
  template void spline_to_model_rotate<dfloat>(std::vector<dfloat> &model, const std::vector<g_vec3<dfloat>> &spline, g_vec3<dfloat> axis, int rotations, bool only_pos);

  template<typename float_type>
  int step_by_spline(const std::vector<g_vec3<float_type>> &spline, g_vec3<float_type> vec)
  {
    int step = 0;
    for (step = 0; step < spline.size() / 2; ++step)
    {
      if (vec.y < spline[step].y)
      {
        break;
      }
    }
    return step;
  }

  template<typename float_type>
  int back_step_by_spline(const std::vector<g_vec3<float_type>> &spline, g_vec3<float_type> vec)
  {
    int step = spline.size() - 1;
    for (step = spline.size() - 1; step >= spline.size() / 2; --step)
    {
      if (vec.y < spline[step].y)
      {
        break;
      }
    }
    return step;
  }

  template<typename float_type>
  float_type get_x_by_spline_with_thick(const std::vector<g_vec3<float_type>> &spline, g_vec3<float_type> vec)
  {
    int step = step_by_spline(spline, vec);
    int step2 = back_step_by_spline(spline, vec);
    if (step != 0 && step < spline.size() / 2 && step2 != spline.size() - 1 && step2 >= spline.size() / 2)
    {
      float_type y1 = (spline[step].y - vec.y) / (spline[step].y - spline[step - 1].y);
      float_type y2 = (vec.y - spline[step - 1].y) / (spline[step].y - spline[step - 1].y);
      float_type x1 = y1 * spline[step - 1].x + y2 * spline[step].x;
      y1 = (spline[step2].y - vec.y) / (spline[step2].y - spline[step2 + 1].y);
      y2 = (vec.y - spline[step2 + 1].y) / (spline[step2].y - spline[step2 + 1].y);
      float_type x2 = y1 * spline[step2 + 1].x + y2 * spline[step2].x;
      return -(x1 + x2) / 2;
    }
    return vec.x;
  }

  template<typename float_type>
  void spline_to_model_part_rotate_plus_shift(std::vector<float_type> &model, const std::vector<g_vec3<float_type>> &spline, const std::vector<g_vec3<float_type>> &cup, g_vec3<float_type> axis, float_type beg_angle, float_type part,
                                              int rotations, g_vec3<float_type> shift, g_vec3<float_type> radius_vec, std::vector<float_type> &radiuses, std::vector<float_type> &thickness, bool only_pos)
  {
    g_mat43<float_type> rot_mat = ident<float_type>();
    g_mat43<float_type> first_rot_mat = ident<float_type>();
    float_type angle = 2*part*PI/(rotations);
    rot_mat = rotate(rot_mat, axis, angle);
    float_type ba = beg_angle + 1e-6;
    first_rot_mat = rotate(first_rot_mat, axis, ba);
    int sp_sz = spline.size();
    int prev_size = model.size()/FLOAT_PER_VERTEX;
    model.reserve(prev_size + FLOAT_PER_VERTEX*3*2*(sp_sz-1));
    std::vector<g_vec3<float_type>> verts = spline;
    std::vector<g_vec3<float_type>> prev_verts = spline;
    std::vector<g_vec3<float_type>> radius_vectors(rotations+1, radius_vec);
    std::vector<g_vec3<float_type>> spline_with_thick = cup;
    if constexpr (!std::is_same<float_type, dfloat>::value)
    {
      spline_with_thick = spline_to_closed_curve_thickness<float_type>(cup, 0.025, 1, 0);
    }

    float_type full_len = 1e-9;
    for (int i=1;i<sp_sz;i++)
    {
      full_len += length(verts[i] - verts[i-1]);
      verts[i - 1] = mulp(first_rot_mat, verts[i - 1]);
      prev_verts[i - 1] = mulp(first_rot_mat, prev_verts[i - 1]);
    }
    int off = (verts.size() - 2) * 6;

    verts[sp_sz - 1] = mulp(first_rot_mat, verts[sp_sz - 1]);
    prev_verts[sp_sz - 1] = mulp(first_rot_mat, prev_verts[sp_sz - 1]);
    radius_vec = mulv(first_rot_mat, radius_vec);
    g_vec3<float_type> prev_radius_vec = radius_vec;

    radius_vectors[0] = radius_vec;
    for (int i=1;i<=rotations;i++)
      radius_vectors[i] = mulv(rot_mat, radius_vectors[i-1]);
    for (int i=0;i<=rotations;i++)
      radius_vectors[i] *= radiuses[i];

    std::vector<g_vec3<float_type>> positions;//size 4*(sp_sz-1)*rotations
    positions.reserve(4*(sp_sz-1)*rotations);
    std::vector<g_vec3<float_type>> normals;//size (sp_sz-1)*rotations
    normals.reserve((sp_sz-1)*rotations);
    std::vector<g_vec2<float_type>> tc;//size 4*(sp_sz-1)*rotations
    tc.reserve(4*(sp_sz-1)*rotations);
    std::vector<int> indices;//size (sp_sz-1)*rotations
    indices.reserve((sp_sz-1)*rotations);

    g_vec3<float_type> vec0 = thickness[0]*verts[0] + radiuses[0]*radius_vec+shift;
    g_vec3<float_type> norm = vec0;
    if constexpr (!std::is_same<float_type, dfloat>::value)
    {
      vec0.x = get_x_by_spline_with_thick(spline_with_thick, vec0);
    }
    for (int i = 2; i < verts.size(); ++i)
    {
      g_vec3<float_type> vec1 = thickness[0]*verts[i-1] + radiuses[0]*radius_vec+shift;
      if constexpr (!std::is_same<float_type, dfloat>::value)
      {
        vec1.x = get_x_by_spline_with_thick(spline_with_thick, vec1);
      }
      if (i == 2)
      {
        norm = normalize(cross(vec0, vec1));
      }
      g_vec3<float_type> vec2 = thickness[0]*verts[i] + radiuses[0]*radius_vec+shift;
      if constexpr (!std::is_same<float_type, dfloat>::value)
      {
        vec2.x = get_x_by_spline_with_thick(spline_with_thick, vec2);
      }
      float_type m = 0.001f;
      add_vertex<float_type>(model, prev_size+0+(i-2)*3, vec0, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+1+(i-2)*3, vec1, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+2+(i-2)*3, vec2, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+0+(i-2)*3 + off/2, vec0 + m*norm, norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+1+(i-2)*3 + off/2, vec1 + m*norm, norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+2+(i-2)*3 + off/2, vec2 + m*norm, norm, {0, 0}, only_pos);//заглушка
    }

    {
      g_vec3<float_type> r_next = radius_vectors[1];
      g_vec3<float_type> v1 = normalize(cross(r_next - radius_vectors[0], axis));
      g_vec3<float_type> v2 = axis;
      float_type phi = (2*PI)/(sp_sz-1);
      for (int i=0;i<sp_sz;i++)
      {
        prev_verts[i] = v1*cos(i*phi) + v2*sin(i*phi);
      }      
    }

    for (int sector = 1; sector <= rotations; sector++)
    {
      g_vec3<float_type> r_next = (sector == rotations) ? radius_vectors[rotations] : radius_vectors[sector+1];
      g_vec3<float_type> v1 = normalize(cross(r_next - radius_vectors[sector-1], axis));
      //std::cerr<<"a "<<radius_vectors[sector-1].x<<" "<<radius_vectors[sector-1].y<<" "<<radius_vectors[sector-1].z<<"\n";
      //std::cerr<<"b "<<v1.x<<" "<<v1.y<<" "<<v1.z<<"\n";
      g_vec3<float_type> v2 = axis;
      float_type phi = (2*PI)/(sp_sz-1);
      for (int i=0;i<sp_sz;i++)
      {
        verts[i] = v1*cos(i*phi) + v2*sin(i*phi);
        //std::cerr<<"verts[i]"<<verts[i].x<<" "<<verts[i].y<<" "<<verts[i].z<<"\n";
      }
      radius_vec = mulv(rot_mat, radius_vec);
      float_type prev_len = 0;
      float_type new_len = 0;
      float_type rv_mult = radiuses[sector];
      float_type prv_mult = radiuses[sector-1];
      float_type thick_mult = thickness[sector];
      float_type pthick_mult = thickness[sector-1];
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + length(verts[i]-verts[i-1]) * thick_mult;
        g_vec3<float_type> p1 = thick_mult*verts[i] + rv_mult*radius_vec + shift;
        g_vec3<float_type> p2 = pthick_mult*prev_verts[i] + prv_mult*prev_radius_vec + shift;
        g_vec3<float_type> p3 = thick_mult*verts[i-1] + rv_mult*radius_vec + shift;
        g_vec3<float_type> p4 = pthick_mult*prev_verts[i-1] + prv_mult*prev_radius_vec + shift;
        if constexpr (!std::is_same<float_type, dfloat>::value)
        {
          if (sector == 1)
          {
            p2.x = get_x_by_spline_with_thick(spline_with_thick, p2);
            p4.x = get_x_by_spline_with_thick(spline_with_thick, p4);
          }
          else if (sector == rotations)
          {
            p1.x = get_x_by_spline_with_thick(spline_with_thick, p1);
            p3.x = get_x_by_spline_with_thick(spline_with_thick, p3);
          }
        }
        g_vec3<float_type> v1 = p2 - p1;
        g_vec3<float_type> v2 = p3 - p1; 
        g_vec3<float_type> n = normalize(cross(v1, v2));
        //std::cerr<<"v1[i]"<<p1.x<<" "<<p1.y<<" "<<p1.z<<"\n";
        //std::cerr<<"v2[i]"<<p3.x<<" "<<p3.y<<" "<<p3.z<<"\n";
        //std::cerr<<"n[i]"<<n.x<<" "<<n.y<<" "<<n.z<<"\n";

        float_type d = dot(n, thick_mult*verts[i]);

        positions.push_back(p1);
        positions.push_back(p2);
        positions.push_back(p3);
        positions.push_back(p4);

        normals.push_back(n);

        tc.push_back(g_vec2<float_type>{((float)(sector))/rotations, 0.75 + 0.25*i/sp_sz});
        tc.push_back(g_vec2<float_type>{((float)(sector - 1))/rotations,  0.75 + 0.25*i/sp_sz});
        tc.push_back(g_vec2<float_type>{((float)(sector))/rotations,  0.75 + 0.25*(i-1)/sp_sz});
        tc.push_back(g_vec2<float_type>{((float)(sector - 1))/rotations,   0.75 + 0.25*(i-1)/sp_sz});

        indices.push_back(prev_size + 6*((sp_sz-1)*(sector-1) + i - 1) + off);

        prev_len = new_len;
      }
      prev_verts = verts;
      prev_radius_vec = radius_vec;
    }

    vec0 = thickness[rotations]*verts[0] + radiuses[rotations]*radius_vec+shift;
    norm = vec0;
    if constexpr (!std::is_same<float_type, dfloat>::value)
    {
      vec0.x = get_x_by_spline_with_thick(spline_with_thick, vec0);
    }
    for (int i = 2; i < verts.size(); ++i)
    {
      g_vec3<float_type> vec1 = thickness[rotations]*verts[i-1] + radiuses[rotations]*radius_vec+shift;
      if constexpr (!std::is_same<float_type, dfloat>::value)
      {
        vec1.x = get_x_by_spline_with_thick(spline_with_thick, vec1);
      }
      if (i == 2)
      {
        norm = normalize(cross(vec0, vec1));
      }
      float_type m = 0.001f;
      g_vec3<float_type> vec2 = thickness[rotations]*verts[i] + radiuses[rotations]*radius_vec+shift;
      if constexpr (!std::is_same<float_type, dfloat>::value)
      {
        vec2.x = get_x_by_spline_with_thick(spline_with_thick, vec2);
      }
      add_vertex<float_type>(model, prev_size+0+(i-2)*3 + off/2, vec0, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+1+(i-2)*3 + off/2, vec1, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+2+(i-2)*3 + off/2, vec2, -norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+0+(i-2)*3 + off/2, vec0 + m*norm, norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+1+(i-2)*3 + off/2, vec1 + m*norm, norm, {0, 0}, only_pos); 
      add_vertex<float_type>(model, prev_size+2+(i-2)*3 + off/2, vec2 + m*norm, norm, {0, 0}, only_pos);//заглушка 2
    }

    for (int i=0;i<(sp_sz-1)*rotations;i++)
    {
      int sp_n = i % (sp_sz-1) + 1;
      int rot_n = i / (sp_sz-1) + 1;

      g_vec3<float_type> n = normals[i];
      g_vec3<float_type> n01 = (sp_n == sp_sz-1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n_11 = (sp_n == sp_sz-1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n11 = (sp_n == sp_sz-1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n];
      g_vec3<float_type> n_10 = normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 1];
      g_vec3<float_type> n10 = normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 1];
      g_vec3<float_type> n0_1 = (sp_n == 1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n - 2];
      g_vec3<float_type> n_1_1 = (sp_n == 1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 2];
      g_vec3<float_type> n1_1 = (sp_n == 1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 2];

      g_vec3<float_type> n0 = normalize(n + n01 + n11 + n10);
      g_vec3<float_type> n1 = normalize(n + n01 + n_11 + n_10);
      g_vec3<float_type> n2 = normalize(n + n10 + n1_1 + n0_1);
      g_vec3<float_type> n3 = normalize(n + n0_1 + n_1_1 + n_10);

      if constexpr (std::is_same<float_type, float>::value)
      {
        bool wrong_norm = false;
        std::vector<g_vec3<float>> norms = {n,n0_1,n_1_1,n_10,n01,n11, n10, n_11, n1_1};
        for (int i0=0;i0<8;i0++)
        {
          for (int j=0;j<8;j++)
          {
            if (i0 != j)
            {
              float d = dot(norms[i0], norms[j]);
              if (d < 0)
                wrong_norm = true;
            }
          }
        }
        if (wrong_norm)
        {
          n0 = n;
          n1 = n;
          n2 = n;
          n3 = n;
        }
      }

      add_vertex<float_type>(model, indices[i]+0, positions[4*i+0], n0, tc[4*i+0], only_pos); 
      add_vertex<float_type>(model, indices[i]+1, positions[4*i+1], n1, tc[4*i+1], only_pos); 
      add_vertex<float_type>(model, indices[i]+2, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex<float_type>(model, indices[i]+3, positions[4*i+3], n3, tc[4*i+3], only_pos); 
      add_vertex<float_type>(model, indices[i]+4, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex<float_type>(model, indices[i]+5, positions[4*i+1], n1, tc[4*i+1], only_pos); //
    }
  }
  template void spline_to_model_part_rotate_plus_shift<dfloat>(std::vector<dfloat> &model, const std::vector<g_vec3<dfloat>> &spline, const std::vector<g_vec3<dfloat>> &cup, g_vec3<dfloat> axis, dfloat beg_angle, dfloat part,
                                              int rotations, g_vec3<dfloat> shift, g_vec3<dfloat> radius_vec, std::vector<dfloat> &radiuses, std::vector<dfloat> &thickness, bool only_pos);

  template<typename float_type>
  PartOffsets create_cup(const std::vector<float_type> &params, std::vector<float_type> &vert, ModelQuality quality)
  {
    int q_pow = (int)pow(2, (int)quality.quality_level);
    int spline_real_start = 1;
    std::vector<g_vec3<float_type>> spline = create_spline(params, 9, 1, 0, true, q_pow);
    if (q_pow > 0)
      spline = spline_make_smoother(spline, 2, spline_real_start, -1, 1, 0);
    g_mat43<float_type> sc = scale(ident<float_type>(), g_vec3<float_type>{0.09,0.9,0.09});
    transform(spline, sc);
    if (params[10] > 0.5)
    {
      int handle_param_idx = 11;
      int radiuses_cnt = 20;
      std::vector<g_vec3<float_type>> spline1 = create_spline_for_handle<float_type>(0, 1);
      float_type thick = 0.02;
      int radius_samples = radiuses_cnt - 1;
      std::vector<float_type> radiuses(radiuses_cnt,0);
      std::vector<float_type> thickness(radiuses_cnt,0);
      for (int i=0;i<radiuses_cnt;i++)
      {
        radiuses[i] = params[handle_param_idx+2 + i];
        thickness[i] = params[handle_param_idx]*params[handle_param_idx+2+radiuses_cnt + i];
      }
      if (radius_samples != radiuses_cnt*2 - 1)
      {
        int k = (radius_samples + 1)/radiuses_cnt;
        std::vector<g_vec3<float_type>> rad_spline = create_spline(radiuses, radiuses_cnt, 0, 1, false, 1);
        rad_spline = spline_make_smoother(rad_spline, k, 0, -1, 0, 1);
        radius_samples = rad_spline.size() - 1;
        radiuses = std::vector<float_type>(rad_spline.size());
        for (int i=0;i<rad_spline.size();i++)
          radiuses[i] = rad_spline[i][1];
      }
      float_type center = params[handle_param_idx+1];
      float_type start_pos = center + radiuses[0];
      float_type end_pos = center - radiuses.back();
      spline1 = spline_rotation(spline1, g_vec3<float_type>{1, 0, 0}, 6*q_pow);
      float_type sin_p;
      if constexpr (std::is_same<float_type, dfloat>::value)
      {
        sin_p = d_min(sin_by_points(spline, start_pos, center, thick, spline_real_start), 0.98);
      }
      else
      {
        sin_p = MIN(sin_by_points(spline, start_pos, center, thick, spline_real_start), 0.98);
      }
      radiuses[0] = dist(x_for_spline_y(spline, center, spline_real_start), center, x_for_spline_y(spline, start_pos, spline_real_start), start_pos);
      radiuses.back() = dist(x_for_spline_y(spline, center, spline_real_start), center, x_for_spline_y(spline, end_pos, spline_real_start), end_pos);
      g_vec3<float_type> radius_vec = g_vec3<float_type>{0, 1, 0};
      spline_to_model_part_rotate_plus_shift<float_type>(vert, spline1, spline, g_vec3<float_type>{0, 0, 1}, asin(sin_p), 0.5, radius_samples, 
                                             shift_by_points(spline, center, thick, 0, 1), 
                                             radius_vec, radiuses, thickness,
                                             quality.create_only_position);
    }
    spline = spline_to_closed_curve_thickness<float_type>(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, spline, g_vec3<float_type>{0,1,0}, 12*q_pow, quality.create_only_position);
    g_mat43<float_type> sc2 = scale(ident<float_type>(), g_vec3<float_type>{1,params[9],1});
    sc2 = translate(sc2, g_vec3<float_type>{0, -0.5,0});
    transform(vert, sc2);

    return dgen::simple_mesh();
  }
  template PartOffsets create_cup<dfloat>(const std::vector<dfloat> &params, std::vector<dfloat> &vert, ModelQuality quality);
  template PartOffsets create_cup<float>(const std::vector<float> &params, std::vector<float> &vert, ModelQuality quality);

  dfloat parameters_cup_reg(const std::vector<dfloat> &params)
  {
    int spline_offsets_cnt = 9;
    dfloat res = 0;

    for (int i = 1; i < spline_offsets_cnt; i++)
    {
      res += 5*d_max(params[i-1] - params[i], 0);
    }
    res = res + 5*d_max(params[spline_offsets_cnt - 1]/(params[0] + 0.001) - 3, 0);
    int handle_param_idx = 11;
    if (params[handle_param_idx - 1] > 0.5)
    {
      int radiuses_cnt = 20;
      if (params[handle_param_idx + 2] <= params[handle_param_idx+radiuses_cnt + 2])
      {
        res += 5;
      }
      if (params[handle_param_idx + 3] <= params[handle_param_idx+radiuses_cnt + 3])
      {
        res += 5;
      }
      for (int i=2;i<radiuses_cnt-1;i++)
      {
        if (params[handle_param_idx + 2 + i] <= params[handle_param_idx+radiuses_cnt + 2 + i])
        {
          res += 5;
        }
        res += 10*d_max(abs(params[handle_param_idx + 2 + i] - params[handle_param_idx + 2 + i - 1]) - 0.1, 0);
        res += 10*d_max(abs(params[handle_param_idx+radiuses_cnt + 2 + i] - params[handle_param_idx+radiuses_cnt + 2 + i - 1]) - 0.5, 0);
      }
    }
    res = d_max(res, 0);
    return res;
  }
};