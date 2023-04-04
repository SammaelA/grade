#include "dishes_generator.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  dfloat triangle_func(dfloat y, int n, int i)
  {
    return d_max(1 - abs(n * y - i), 0);
  }

  dfloat x_for_spline_y(const std::vector<dvec3> &in_spline, dfloat y, int start_id)
  {
    dfloat sum = 0;
    if (y < in_spline[0][1])
      return in_spline[0][0];
    if (y > in_spline.back()[1])
      return in_spline.back()[0];

    for (int i = start_id + 1; i < in_spline.size(); ++i)
    {
      if (y >= in_spline[i-1][1] && y <= in_spline[i][1])
      {
        dfloat d = in_spline[i][1] - in_spline[i-1][1];
        dfloat d1 = ((y - in_spline[i-1][1])/d);
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

  dfloat dist(dfloat x1, dfloat y1, dfloat x2, dfloat y2)
  {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  }

  dfloat rad_by_points(const std::vector<dvec3> &in_spline, dfloat y1, dfloat y2, int start_id)
  {
    return dist(x_for_spline_y(in_spline, y1, start_id), y1, x_for_spline_y(in_spline, y2, start_id), y2) / 2.0;
  }

  dvec3 shift_by_points(const std::vector<dvec3> &in_spline, dfloat center_y, dfloat thick, int x, int y)
  {
    dvec3 shift{0, 0, 0};
    shift[y] = center_y;
    shift[x] = -x_for_spline_y(in_spline, center_y, 0) + thick;
    return shift;
  }

  dfloat sin_by_points(const std::vector<dvec3> &in_spline, dfloat y1, dfloat y2, dfloat thick, int start_id)
  {
    return (x_for_spline_y(in_spline, y1, start_id) - x_for_spline_y(in_spline, y2, start_id) + thick / 2.0) / (2.0 * rad_by_points(in_spline, y1, y2, start_id));
  }

  std::vector<dvec3> create_spline(const std::vector<dfloat> &params, int idx, int axis_x, int axis_y, bool from_zero, int detail_q)
  {
    //if you have n params and axis_x = 0, axis_y = 2 then it will create n points (vec3) with formula point[i] = (i/(n-1), 0, params[i])
    std::vector<dvec3> spline;
    spline.reserve((params.size() + from_zero));
    if (from_zero)
    {
      dfloat smooth_y = d_min(0.5, 0.9*params[0]);
      dfloat smooth_x = 0.05;
      //center of cup bottom
      dvec3 vec{0,0,0};
      vec[axis_x] = -smooth_x -1e-4;
      spline.push_back(vec);

      //smoothness on the bottom/side transition
      int steps = 2 + 2*detail_q;
      for (int i=0;i<steps;i++)
      {
        dvec3 vec{0,0,0};
        vec[axis_x] = -smooth_x*cos((PI*i)/(2*steps));
        vec[axis_y] = params[0] - smooth_y*(1 - sin((PI*i)/(2*steps)));
        spline.push_back(vec);
      }
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
    dmat43 rot_mat = ident<dfloat>();
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

  std::vector<dvec3> spline_to_closed_curve_thickness(const std::vector<dvec3> &in_spline, dfloat thickness, int axis_x, int axis_y)
  {
    std::vector<dvec3> spline = in_spline;
    std::vector<dvec3> offset_dirs;
    dvec3 x = dvec3{0,0,0};
    x[axis_x] = 1;
    dvec3 y = dvec3{0,0,0};
    y[axis_y] = -1;
    offset_dirs.push_back(in_spline[0] + thickness*x);

    dvec3 l0 = in_spline[1] - in_spline[0];
    dvec3 d_prev = dvec3{0,0,0};
    d_prev[axis_x] = l0[axis_y];
    d_prev[axis_y] = -l0[axis_x];
    for (int i=1;i<in_spline.size()-1;i++)
    {
      dvec3 l = in_spline[i+1] - in_spline[i];
      dvec3 d = dvec3{0,0,0};
      d[axis_x] = l[axis_y];
      d[axis_y] = -l[axis_x];
      offset_dirs.push_back(in_spline[i] + thickness*normalize(d_prev + d)); //pos + thickness*dir
      d_prev = d;
    }
    offset_dirs.push_back(in_spline[in_spline.size()-1] + thickness*y);

    //smooth transition on top
    dfloat d = spline.back()[axis_y] - offset_dirs.back()[axis_y];
    {
      dvec3 v = spline.back();
      v[axis_x] += 0.25*d;
      v[axis_y] -= 0.25*d;
      spline.push_back(v);
    }
    {
      dvec3 v = spline.back();
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

  dfloat lagrange_3_dots(dfloat x, dfloat x1, dfloat y1, dfloat x2, dfloat y2, dfloat x3, dfloat y3)
  {
    dfloat y =  y1*((x - x2)/(x1 - x2))*((x - x3)/(x1 - x3)) +
           y2*((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3)) +
           y3*((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2));
    //std::cerr<<x<<"--"<<y<<"--"<<x1<<"--"<<x2<<"--"<<x3<<"--"<<((x - x2)/(x1 - x2))*((x - x3)/(x1 - x2))<<"--"<<((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3))<<"--"<<((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2))<<"\n";
    return y;
  }

  std::vector<dvec3> spline_make_smoother(const std::vector<dvec3> &in_spline, int detail_q, int index_from, int index_to, int axis_x, int axis_y)
  {
    index_from = index_from < 0 ? 0 : index_from;
    index_to = index_to < 0 ? in_spline.size() : index_to;
    if (detail_q <= 1)
      return in_spline;
    if (detail_q == 2 && (in_spline.size() == index_to))
    {
      std::vector<dvec3> sp = std::vector<dvec3>(index_from + 2*(index_to - index_from) - 1, dvec3{0,0,0});
      for (int i = 0; i < index_from; i++)
        sp[i] = in_spline[i];
      for (int i = index_from; i < index_to; i++)
        sp[index_from + 2*(i - index_from)] = in_spline[i];
      for (int i = index_from+1; i < index_to-1; i++)
      {
        dfloat new_x1 = (in_spline[i-1][axis_x] + in_spline[i][axis_x])/2;
        dfloat new_x2 = (in_spline[i][axis_x] + in_spline[i+1][axis_x])/2;
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
    dmat43 rot_mat = ident<dfloat>();
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
      full_len += length(verts[i]-verts[i-1]);
    }

    std::vector<dvec3> positions;//size 4*(sp_sz-1)*rotations
    positions.reserve(4*(sp_sz-1)*rotations);
    std::vector<dvec3> normals;//size (sp_sz-1)*rotations
    normals.reserve((sp_sz-1)*rotations);
    std::vector<dvec2> tc;//size 4*(sp_sz-1)*rotations
    tc.reserve(4*(sp_sz-1)*rotations);
    std::vector<int> indices;//size (sp_sz-1)*rotations
    indices.reserve((sp_sz-1)*rotations);

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
        new_len = prev_len + length(verts[i] - verts[i-1]);
        dvec3 v1 = prev_verts[i-1] - verts[i];
        dvec3 v2 = verts[i-1] - prev_verts[i]; 
        dvec3 n = normalize(cross(v1, v2));

        positions.push_back(verts[i]);
        positions.push_back(prev_verts[i]);
        positions.push_back(verts[i-1]);
        positions.push_back(prev_verts[i-1]);

        normals.push_back(n);

        tc.push_back(dvec2{((float)(sector))/rotations, 0.75*new_len/full_len});
        tc.push_back(dvec2{((float)(sector - 1))/rotations, 0.75*new_len/full_len});
        tc.push_back(dvec2{((float)(sector))/rotations, 0.75*prev_len/full_len});
        tc.push_back(dvec2{((float)(sector - 1))/rotations,  0.75*prev_len/full_len});

        indices.push_back(prev_size + 6*((sp_sz-1)*(sector-1) + i - 1));
        
        prev_len = new_len;
      }
      prev_verts = verts;
    }

    for (int i=0;i<(sp_sz-1)*rotations;i++)
    {
      int sp_n = i % (sp_sz-1) + 1;
      int rot_n = i / (sp_sz-1) + 1;

      dvec3 n = normals[i];
      dvec3 n01 = (sp_n == sp_sz-1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n];
      dvec3 n_11 = (sp_n == sp_sz-1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n];
      dvec3 n11 = (sp_n == sp_sz-1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n];
      dvec3 n_10 = normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 1];
      dvec3 n10 = normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 1];
      dvec3 n0_1 = (sp_n == 1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n - 2];
      dvec3 n_1_1 = (sp_n == 1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 2];
      dvec3 n1_1 = (sp_n == 1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 2];

      dvec3 n0 = normalize(n + n01 + n11 + n10);
      dvec3 n1 = normalize(n + n01 + n_11 + n_10);
      dvec3 n2 = normalize(n + n10 + n1_1 + n0_1);
      dvec3 n3 = normalize(n + n0_1 + n_1_1 + n_10);
      add_vertex(model, indices[i]+0, positions[4*i+0], n0, tc[4*i+0], only_pos); 
      add_vertex(model, indices[i]+1, positions[4*i+1], n1, tc[4*i+1], only_pos); 
      add_vertex(model, indices[i]+2, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex(model, indices[i]+3, positions[4*i+3], n3, tc[4*i+3], only_pos); 
      add_vertex(model, indices[i]+4, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex(model, indices[i]+5, positions[4*i+1], n1, tc[4*i+1], only_pos); 
    }
  }

  void spline_to_model_part_rotate_plus_shift(std::vector<dfloat> &model, const std::vector<dvec3> &spline, dvec3 axis, dfloat beg_angle, dfloat part,
                                              int rotations, dvec3 shift, dvec3 radius_vec, std::vector<dfloat> &radiuses, std::vector<dfloat> &thickness, bool only_pos)
  {
    dmat43 rot_mat = ident<dfloat>();
    dmat43 first_rot_mat = ident<dfloat>();
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
      full_len += length(verts[i] - verts[i-1]);
      verts[i - 1] = mulp(first_rot_mat, verts[i - 1]);
      prev_verts[i - 1] = mulp(first_rot_mat, prev_verts[i - 1]);
    }

    verts[sp_sz - 1] = mulp(first_rot_mat, verts[sp_sz - 1]);
    prev_verts[sp_sz - 1] = mulp(first_rot_mat, prev_verts[sp_sz - 1]);
    radius_vec = mulv(first_rot_mat, radius_vec);
    dvec3 prev_radius_vec = radius_vec;

    std::vector<dvec3> positions;//size 4*(sp_sz-1)*rotations
    positions.reserve(4*(sp_sz-1)*rotations);
    std::vector<dvec3> normals;//size (sp_sz-1)*rotations
    normals.reserve((sp_sz-1)*rotations);
    std::vector<dvec2> tc;//size 4*(sp_sz-1)*rotations
    tc.reserve(4*(sp_sz-1)*rotations);
    std::vector<int> indices;//size (sp_sz-1)*rotations
    indices.reserve((sp_sz-1)*rotations);

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
      dfloat thick_mult = thickness[sector];
      dfloat pthick_mult = thickness[sector-1];
      for (int i=1;i<sp_sz;i++)
      {
        new_len = prev_len + length(verts[i]-verts[i-1]) * thick_mult;
        dvec3 p1 = thick_mult*verts[i] + rv_mult*radius_vec + shift;
        dvec3 p2 = pthick_mult*prev_verts[i] + prv_mult*prev_radius_vec + shift;
        dvec3 p3 = thick_mult*verts[i-1] + rv_mult*radius_vec + shift;
        dvec3 p4 = pthick_mult*prev_verts[i-1] + prv_mult*prev_radius_vec + shift;
        dvec3 v1 = p2 - p1;
        dvec3 v2 = p3 - p1; 
        dvec3 n = normalize(cross(v1, v2));

        positions.push_back(p1);
        positions.push_back(p2);
        positions.push_back(p3);
        positions.push_back(p4);

        normals.push_back(n);

        tc.push_back(dvec2{((float)(sector))/rotations, 0.75 + 0.25*new_len/full_len});
        tc.push_back(dvec2{((float)(sector - 1))/rotations,  0.75 + 0.25*new_len/full_len});
        tc.push_back(dvec2{((float)(sector))/rotations,  0.75 + 0.25*prev_len/full_len});
        tc.push_back(dvec2{((float)(sector - 1))/rotations,   0.75 + 0.25*prev_len/full_len});

        indices.push_back(prev_size + 6*((sp_sz-1)*(sector-1) + i - 1));

        prev_len = new_len;
      }
      prev_verts = verts;
      prev_radius_vec = radius_vec;
    }

    for (int i=0;i<(sp_sz-1)*rotations;i++)
    {
      int sp_n = i % (sp_sz-1) + 1;
      int rot_n = i / (sp_sz-1) + 1;

      dvec3 n = normals[i];
      dvec3 n01 = (sp_n == sp_sz-1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n];
      dvec3 n_11 = (sp_n == sp_sz-1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n];
      dvec3 n11 = (sp_n == sp_sz-1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n];
      dvec3 n_10 = normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 1];
      dvec3 n10 = normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 1];
      dvec3 n0_1 = (sp_n == 1) ? n : normals[(rot_n-1)*(sp_sz-1) + sp_n - 2];
      dvec3 n_1_1 = (sp_n == 1) ? n : normals[((rot_n-2 + rotations) % rotations)*(sp_sz-1) + sp_n - 2];
      dvec3 n1_1 = (sp_n == 1) ? n : normals[(rot_n % rotations)*(sp_sz-1) + sp_n - 2];

      dvec3 n0 = normalize(n + n01 + n11 + n10);
      dvec3 n1 = normalize(n + n01 + n_11 + n_10);
      dvec3 n2 = normalize(n + n10 + n1_1 + n0_1);
      dvec3 n3 = normalize(n + n0_1 + n_1_1 + n_10);
      add_vertex(model, indices[i]+0, positions[4*i+0], n0, tc[4*i+0], only_pos); 
      add_vertex(model, indices[i]+1, positions[4*i+1], n1, tc[4*i+1], only_pos); 
      add_vertex(model, indices[i]+2, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex(model, indices[i]+3, positions[4*i+3], n3, tc[4*i+3], only_pos); 
      add_vertex(model, indices[i]+4, positions[4*i+2], n2, tc[4*i+2], only_pos); 
      add_vertex(model, indices[i]+5, positions[4*i+1], n1, tc[4*i+1], only_pos); 
    }
  }

  PartOffsets create_cup(const std::vector<dfloat> &params, std::vector<dfloat> &vert, ModelQuality quality)
  {
    int q_pow = (int)pow(2, (int)quality.quality_level);
    int spline_real_start = 1;
    std::vector<dvec3> spline = create_spline(params, 9, 1, 0, true, q_pow);
    if (q_pow > 0)
      spline = spline_make_smoother(spline, 2, spline_real_start, -1, 1, 0);
    dmat43 sc = scale(ident<dfloat>(), dvec3{0.09,0.9,0.09});
    transform(spline, sc);
    if (params[10] > 0.5)
    {
      int handle_param_idx = 11;
      int radiuses_cnt = 20;
      std::vector<dvec3> spline1 = create_spline_for_handle(params, handle_param_idx, 0, 1);
      dfloat thick = 0.02;
      int radius_samples = radiuses_cnt - 1;
      std::vector<dfloat> radiuses(radiuses_cnt,0);
      std::vector<dfloat> thickness(radiuses_cnt,0);
      for (int i=0;i<radiuses_cnt;i++)
      {
        radiuses[i] = params[handle_param_idx+2 + i];
        thickness[i] = params[handle_param_idx+2+radiuses_cnt + i];
      }
      if (radius_samples != radiuses_cnt*2 - 1)
      {
        int k = (radius_samples + 1)/radiuses_cnt;
        std::vector<dvec3> rad_spline = create_spline(radiuses, radiuses_cnt, 0, 1, false, 1);
        rad_spline = spline_make_smoother(rad_spline, k, 0, -1, 0, 1);
        radius_samples = rad_spline.size() - 1;
        radiuses = std::vector<dfloat>(rad_spline.size());
        for (int i=0;i<rad_spline.size();i++)
          radiuses[i] = rad_spline[i][1];
      }
      dfloat center = params[handle_param_idx+1];
      dfloat start_pos = center + radiuses[0];
      dfloat end_pos = center - radiuses.back();
      spline1 = spline_rotation(spline1, dvec3{1, 0, 0}, 6*q_pow);
      dfloat sin_p = d_min(sin_by_points(spline, start_pos, center, thick, spline_real_start), 0.98);
      radiuses[0] = dist(x_for_spline_y(spline, center, spline_real_start), center, x_for_spline_y(spline, start_pos, spline_real_start), start_pos);
      radiuses.back() = dist(x_for_spline_y(spline, center, spline_real_start), center, x_for_spline_y(spline, end_pos, spline_real_start), end_pos);
      dvec3 radius_vec = dvec3{0, 1, 0};
      spline_to_model_part_rotate_plus_shift(vert, spline1, dvec3{0, 0, 1}, asin(sin_p), 0.5, radius_samples, 
                                             shift_by_points(spline, center, thick, 0, 1), 
                                             radius_vec, radiuses, thickness,
                                             quality.create_only_position);
    }
    spline = spline_to_closed_curve_thickness(spline, 0.025, 1, 0);
    spline_to_model_rotate(vert, spline, dvec3{0,1,0}, 12*q_pow, quality.create_only_position);
    dmat43 sc2 = scale(ident<dfloat>(), dvec3{1,params[9],1});
    sc2 = translate(sc2, dvec3{0, -0.5,0});
    transform(vert, sc2);

    return dgen::simple_mesh();
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

    int handle_param_idx = 11;
    if (params[handle_param_idx - 1] > 0.5)
    {
      int radiuses_cnt = 20;
      for (int i=1;i<radiuses_cnt;i++)
        res += 10*d_max(abs(params[handle_param_idx + 2 + i] - params[handle_param_idx + 2 + i - 1]) - 0.1, 0);
    }

    res = d_max(res, 0);
    return res;
  }
};