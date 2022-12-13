#include "dishes_generator.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
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

  dvec3 shift_by_points(const std::vector<dvec3> &in_spline, dfloat center_y, dfloat thick, int x, int y)
  {
    dvec3 shift{0, 0, 0};
    shift[y] = center_y;
    shift[x] = -x_for_spline_y(in_spline, center_y, 0) + thick;
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

  void spline_to_model_rotate(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, int rotations, bool only_pos)
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
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1), verts[i], n, dvec2{((float)(sector))/rotations, 0.75*new_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+1, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, 0.75*new_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+2, verts[i-1], n, dvec2{((float)(sector))/rotations, 0.75*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+3, prev_verts[i-1], n, dvec2{((float)(sector - 1))/rotations,  0.75*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+4, verts[i-1], n, dvec2{((float)(sector))/rotations,  0.75*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+5, prev_verts[i], n, dvec2{((float)(sector - 1))/rotations, 0.75*new_len/full_len}, only_pos); 
        prev_len = new_len;
      }
      prev_verts = verts;
    }
  }

  void spline_to_model_part_rotate_plus_shift(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, dfloat beg_angle, dfloat part,
                                              int rotations, dvec3 shift, dvec3 radius_vec, std::vector<dfloat> &radiuses, bool only_pos)
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
    radius_vec = mulv(first_rot_mat, radius_vec);
    dvec3 prev_radius_vec = radius_vec;

    for (int sector = 1; sector <= rotations; sector++)
    {
      for (int i=0;i<sp_sz;i++)
      {
        verts[i] = mulp(rot_mat, verts[i]);
      }
      radius_vec = mulv(rot_mat, radius_vec);
      dfloat prev_len = 0;
      dfloat new_len = 0;
      dfloat rv_mult = radiuses[sector];
      dfloat prv_mult = radiuses[sector-1];
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + len(sub(verts[i],verts[i-1]));
        dvec3 p1 = add(add(verts[i],   mul(rv_mult, radius_vec)), shift);
        dvec3 p2 = add(add(prev_verts[i],   mul(prv_mult, prev_radius_vec)), shift);
        dvec3 p3 = add(add(verts[i-1], mul(rv_mult, radius_vec)), shift);
        dvec3 p4 = add(add(prev_verts[i-1], mul(prv_mult, prev_radius_vec)), shift);
        dvec3 v1 = sub(p2, p1);
        dvec3 v2 = sub(p3, p1); 
        dvec3 n = normalize(cross(v1, v2));

        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1),   p1, n, dvec2{((float)(sector))/rotations, 0.75+0.25*new_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+1, p2, n, dvec2{((float)(sector - 1))/rotations, 0.75+0.25*new_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+2, p3, n, dvec2{((float)(sector))/rotations, 0.75+0.25*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+3, p4, n, dvec2{((float)(sector - 1))/rotations,  0.75+0.25*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+4, p3, n, dvec2{((float)(sector))/rotations,  0.75+0.25*prev_len/full_len}, only_pos); 
        add_vertex(model, prev_size + 6*((sp_sz-1)*(sector-1) + i - 1)+5, p2, n, dvec2{((float)(sector - 1))/rotations, 0.75+0.25*new_len/full_len}, only_pos); 
        prev_len = new_len;
      }
      prev_verts = verts;
      prev_radius_vec = radius_vec;
    }
  }

  void create_cup(const std::vector<dfloat> &params, std::vector<dfloat> &vert, ModelQuality quality)
  {
    int q_pow = (int)pow(2, (int)quality.quality_level);
    std::vector<dvec3> spline = create_spline(params, 9, 1, 0, true);
    //if (q_pow > 1)
    //  spline = spline_make_smoother(spline, q_pow, 1, -1, 1, 0);
    dmat43 sc = scale(ident(), dvec3{0.09,0.9,0.09});
    transform(spline, sc);

    if (params[10] > 0.5)
    {
      int handle_param_idx = 11;
      int radiuses_cnt = 7;
      std::vector<dvec3> spline1 = create_spline_for_handle(params, handle_param_idx, 0, 1);
      dfloat thick = smoothmin(params[handle_param_idx], 0.02, 8);
      int radius_samples = radiuses_cnt*q_pow - 1;
      std::vector<dfloat> radiuses(radiuses_cnt,0);
      for (int i=0;i<radiuses_cnt;i++)
        radiuses[i] = params[handle_param_idx+2 + i] + params[handle_param_idx];
      if (radius_samples != radiuses_cnt - 1)
      {
        int k = (radius_samples + 1)/radiuses_cnt;
        std::vector<dvec3> rad_spline = create_spline(radiuses, radiuses_cnt, 0, 1, false);
        //rad_spline = spline_make_smoother(rad_spline, k, 0, -1, 0, 1);
        radius_samples = rad_spline.size() - 1;
        radiuses = std::vector<dfloat>(rad_spline.size());
        for (int i=0;i<rad_spline.size();i++)
          radiuses[i] = rad_spline[i][1];
      }
      dfloat center = params[handle_param_idx+1];
      dfloat start_pos = center + radiuses[0];
      dfloat end_pos = center - radiuses.back();

      spline1 = spline_rotation(spline1, dvec3{1, 0, 0}, 4*q_pow);
      dfloat sin_p = smoothmin(sin_by_points(spline, start_pos, center, thick), 0.98, 8);
      radiuses[0] = dist(x_for_spline_y(spline, center, 0), center, x_for_spline_y(spline, start_pos, 0), start_pos);
      radiuses.back() = dist(x_for_spline_y(spline, center, 0), center, x_for_spline_y(spline, end_pos, 0), end_pos);
      dvec3 radius_vec = dvec3{0, 1, 0};
      spline_to_model_part_rotate_plus_shift(vert, spline1, dvec3{0, 0, 1}, asin(sin_p), 0.5, radius_samples, 
                                             shift_by_points(spline, center, thick, 0, 1), 
                                             radius_vec, radiuses,
                                             quality.create_only_position);
    }
    spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, spline, dvec3{0,1,0}, 12*q_pow, quality.create_only_position);
    dmat43 sc2 = scale(ident(), dvec3{1,params[9],1});
    sc2 = translate(sc2, dvec3{0, -0.5,0});
    transform(vert, sc2);
  }

  dfloat parameters_cup_reg(const std::vector<dfloat> &params)
  {
    int spline_offsets_cnt = 9;
    dfloat res = 0;
    for (int i = 1; i < spline_offsets_cnt; i++)
    {
      res += 5*d_max(params[i-1] - params[i], 0);
    }
    res = res + 5*d_max(params[spline_offsets_cnt - 1]/(params[0] + 0.001) - 3, 0);
    res = d_max(res, 0);
    return res;
  }
};