#include "building_generator.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  std::vector<dfloat> create_wall_spline(const std::vector<dfloat> &params, int idx)
  {
    std::vector<dfloat> spline;
    spline.push_back(0);
    spline.push_back(params[idx]);
    spline.push_back(spline[spline.size() - 1] + params[idx + 3]);
    for (dfloat len = params[idx] + params[idx + 2] + params[idx + 3]; len <= params[idx + 4]; len += params[idx + 3] + params[idx + 1])
    {
      spline.push_back(spline[spline.size() - 1] + params[idx + 1]);
      spline.push_back(spline[spline.size() - 1] + params[idx + 3]);
    }
    spline.push_back(spline[spline.size() - 1] + params[idx + 2]);
    return spline;
  }

  void get_walls_from_splines(std::vector<dfloat> &model, const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, dfloat length_z, dfloat window_depth, int x, int y, int z, float tex_v_shift, bool only_pos)
  {
    int num = model.size();
    dvec3 norm_x = {0, 0, 0};
    dvec3 norm_y = {0, 0, 0};
    norm_x[x] = 1;
    norm_y[y] = 1;
    for (int i = 1; i < sp_x.size(); ++i)
    {
      for (int j = 1; j < sp_y.size(); ++j)
      {
        for (int k = 0; k <= 1; ++k)
        {
          dvec3 depth = {0, 0, 0};

          dvec3 pos_0 = {0, 0, 0};
          dvec3 pos_1 = {0, 0, 0};
          dvec3 pos_2 = {0, 0, 0};
          dvec3 pos_3 = {0, 0, 0};

          dvec2 tex_0 = {0.5 * sp_x[i - 1] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j - 1] / sp_y[sp_y.size() - 1]};
          dvec2 tex_1 = {0.5 * sp_x[i - 1] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j] / sp_y[sp_y.size() - 1]};
          dvec2 tex_2 = {0.5 * sp_x[i] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j - 1] / sp_y[sp_y.size() - 1]};
          dvec2 tex_3 = {0.5 * sp_x[i] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j] / sp_y[sp_y.size() - 1]};

          dvec3 norm_z = {0, 0, 0};
          pos_0[x] = sp_x[i - 1];
          pos_0[y] = sp_y[j - 1];

          pos_1[x] = sp_x[i - 1];
          pos_1[y] = sp_y[j];

          pos_2[x] = sp_x[i];
          pos_2[y] = sp_y[j - 1];

          pos_3[x] = sp_x[i];
          pos_3[y] = sp_y[j];

          norm_z[z] = -1;
          if (k)
          {
            pos_0[z] = length_z;
            pos_1[z] = length_z;
            pos_2[z] = length_z;
            pos_3[z] = length_z;
            tex_0[0] += 0.5;
            tex_1[0] += 0.5;
            tex_2[0] += 0.5;
            tex_3[0] += 0.5;
            norm_z[z] = 1;
          }
          if (i % 2 == 0 && j % 2 == 0)
          {
            depth[z] = window_depth;
            if (k)
            {
              depth[z] = -window_depth;
            }
            dfloat u = 0.2 * 0.5 * (sp_x[i] - sp_x[i - 1]) / sp_x[sp_x.size() - 1];
            dfloat v = 0.2 * 0.5 * (sp_y[j] - sp_y[j - 1]) / sp_y[sp_y.size() - 1];
            
            add_vertex(model, num++, pos_0, norm_x, tex_0 + dvec2{0, v}, only_pos); 
            add_vertex(model, num++, pos_1, norm_x, tex_1 - dvec2{0, v}, only_pos); 
            add_vertex(model, num++, pos_0 + depth, norm_x, tex_0 + dvec2{u, v}, only_pos); 
            add_vertex(model, num++, pos_1 + depth, norm_x, tex_1 + dvec2{u, -v}, only_pos); 
            add_vertex(model, num++, pos_0 + depth, norm_x, tex_0 + dvec2{u, v}, only_pos); 
            add_vertex(model, num++, pos_1, norm_x, tex_1 - dvec2{0, v}, only_pos); 

            add_vertex(model, num++, pos_2, dvec3{0, 0, 0} - norm_x, tex_2 + dvec2{0, v}, only_pos); 
            add_vertex(model, num++, pos_3, dvec3{0, 0, 0} - norm_x, tex_3 - dvec2{0, v}, only_pos); 
            add_vertex(model, num++, pos_2 + depth, dvec3{0, 0, 0} - norm_x, tex_2 + dvec2{-u, v}, only_pos); 
            add_vertex(model, num++, pos_3 + depth, dvec3{0, 0, 0} - norm_x, tex_3 + dvec2{-u, -v}, only_pos); 
            add_vertex(model, num++, pos_2 + depth, dvec3{0, 0, 0} - norm_x, tex_2 + dvec2{-u, v}, only_pos); 
            add_vertex(model, num++, pos_3, dvec3{0, 0, 0} - norm_x, tex_3 - dvec2{0, v}, only_pos);

            add_vertex(model, num++, pos_0, norm_y, tex_0 + dvec2{u, 0}, only_pos); 
            add_vertex(model, num++, pos_2, norm_y, tex_2 - dvec2{u, 0}, only_pos); 
            add_vertex(model, num++, pos_0 + depth, norm_y, tex_0 + dvec2{u, v}, only_pos); 
            add_vertex(model, num++, pos_2 + depth, norm_y, tex_2 + dvec2{-u, v}, only_pos); 
            add_vertex(model, num++, pos_0 + depth, norm_y, tex_0 + dvec2{u, v}, only_pos); 
            add_vertex(model, num++, pos_2, norm_y, tex_2 - dvec2{u, 0}, only_pos); 

            add_vertex(model, num++, pos_1, dvec3{0, 0, 0} - norm_y, tex_1 + dvec2{u, 0}, only_pos); 
            add_vertex(model, num++, pos_3, dvec3{0, 0, 0} - norm_y, tex_3 - dvec2{u, 0}, only_pos); 
            add_vertex(model, num++, pos_1 + depth, dvec3{0, 0, 0} - norm_y, tex_1 + dvec2{u, -v}, only_pos); 
            add_vertex(model, num++, pos_3 + depth, dvec3{0, 0, 0} - norm_y, tex_3 + dvec2{-u, -v}, only_pos); 
            add_vertex(model, num++, pos_1 + depth, dvec3{0, 0, 0} - norm_y, tex_1 + dvec2{u, -v}, only_pos); 
            add_vertex(model, num++, pos_3, dvec3{0, 0, 0} - norm_y, tex_3 - dvec2{u, 0}, only_pos);

            tex_0 += dvec2{u, v};
            tex_1 += dvec2{u, -v};
            tex_2 += dvec2{-u, v};
            tex_3 += dvec2{-u, -v};
          }
          else
          {
            add_vertex(model, num++, pos_0 + depth, norm_z, tex_0, only_pos); 
            add_vertex(model, num++, pos_1 + depth, norm_z, tex_1, only_pos); 
            add_vertex(model, num++, pos_2 + depth, norm_z, tex_2, only_pos); 
            add_vertex(model, num++, pos_3 + depth, norm_z, tex_3, only_pos); 
            add_vertex(model, num++, pos_2 + depth, norm_z, tex_2, only_pos); 
            add_vertex(model, num++, pos_1 + depth, norm_z, tex_1, only_pos); 
          }
        }
      }
    }
  }

  void splines_to_building(std::vector<dfloat> &model, const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, const std::vector<dfloat> &sp_z, dfloat window_depth, bool only_pos)
  {
    get_walls_from_splines(model, sp_x, sp_y, sp_z[sp_z.size() - 1], window_depth, 0, 1, 2, 0, only_pos);
    get_walls_from_splines(model, sp_z, sp_y, sp_x[sp_x.size() - 1], window_depth, 2, 1, 0, 0.35, only_pos);
    int num = model.size();
    add_vertex(model, num++, dvec3{0, 0, 0}, dvec3{0, -1, 0}, dvec2{0, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{0, 0, sp_z[sp_z.size() - 1]}, dvec3{0, -1, 0}, dvec2{0, 1}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], 0, 0}, dvec3{0, -1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], 0, sp_z[sp_z.size() - 1]}, dvec3{0, -1, 0}, dvec2{0.5, 1}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], 0, 0}, dvec3{0, -1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{0, 0, sp_z[sp_z.size() - 1]}, dvec3{0, -1, 0}, dvec2{0, 1}, only_pos); 
    
    add_vertex(model, num++, dvec3{0, sp_y[sp_y.size() - 1], 0}, dvec3{0, 1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{0, sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, dvec3{0, 1, 0}, dvec2{0.5, 1}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], 0}, dvec3{0, 1, 0}, dvec2{1, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, dvec3{0, 1, 0}, dvec2{1, 1}, only_pos); 
    add_vertex(model, num++, dvec3{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], 0}, dvec3{0, 1, 0}, dvec2{1, 0.7}, only_pos); 
    add_vertex(model, num++, dvec3{0, sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, dvec3{0, 1, 0}, dvec2{0.5, 1}, only_pos); 
  }

  void create_building(const std::vector<dfloat> &params, std::vector<dfloat> &vert, ModelQuality quality)
  {
    std::vector<dfloat> spline_x = create_wall_spline(params, 0);
    std::vector<dfloat> spline_y = create_wall_spline(params, 5);
    std::vector<dfloat> spline_z = create_wall_spline(params, 10);
    splines_to_building(vert, spline_x, spline_y, spline_z, params[15], quality.create_only_position);
    dmat43 sc2 = translate(scale(ident(), dvec3{0.02, 0.02, 0.02}), dvec3{-0.5 * params[4], 0, 0});
    transform(vert, sc2);
  }
};