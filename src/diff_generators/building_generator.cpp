#include "building_generator.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  template<typename float_type>
  std::vector<float_type> create_wall_spline(const std::vector<float_type> &params, int idx)
  {
    std::vector<float_type> spline;
    spline.push_back(0);
    spline.push_back(params[idx]);
    spline.push_back(spline[spline.size() - 1] + params[idx + 3]);
    for (float_type len = params[idx] + params[idx + 2] + params[idx + 3]; len <= params[idx + 4]; len += params[idx + 3] + params[idx + 1])
    {
      spline.push_back(spline[spline.size() - 1] + params[idx + 1]);
      spline.push_back(spline[spline.size() - 1] + params[idx + 3]);
    }
    spline.push_back(spline[spline.size() - 1] + params[idx + 2]);
    return spline;
  }
  template std::vector<dfloat> create_wall_spline<dfloat>(const std::vector<dfloat> &params, int idx);

  template<typename float_type>
  std::vector<float_type> create_wall_spline(float_type left_offset, float_type window_gap, float_type right_offset, float_type window_size, 
                                         float_type window_count)
  {
    //std::cerr<<"w count "<<window_count<<"\n";
    std::vector<float_type> spline;
    spline.push_back(0);
    spline.push_back(left_offset);
    if (window_count > 0)
    {
      spline.push_back(spline[spline.size() - 1] + window_size);
      for (float_type i = 0; i < window_count - 1; i += 1)
      {
        spline.push_back(spline[spline.size() - 1] + window_size);
        spline.push_back(spline[spline.size() - 1] + window_size);
      }
    }
    spline.push_back(spline[spline.size() - 1] + right_offset);
    return spline;
  }
  template std::vector<dfloat> create_wall_spline<dfloat>(dfloat left_offset, dfloat window_gap, dfloat right_offset, dfloat window_size, 
                                         dfloat window_count);

  template<typename float_type>
  void get_walls_from_splines(std::vector<float_type> &model, std::vector<float_type> &windows,
                              const std::vector<float_type> &sp_x, const std::vector<float_type> &sp_y, 
                              float_type length_z, float_type window_depth, int x, int y, int z, float tex_v_shift, 
                              bool low_quality, bool only_pos)
  {
    int num = model.size()/FLOAT_PER_VERTEX;
    g_vec3<float_type> norm_x = {0, 0, 0};
    g_vec3<float_type> norm_y = {0, 0, 0};
    norm_x[x] = 1;
    norm_y[y] = 1;
    for (int i = 1; i < sp_x.size(); ++i)
    {
      for (int j = 1; j < sp_y.size(); ++j)
      {
        for (int k = 0; k <= 1; ++k)
        {
          g_vec3<float_type> depth = {0, 0, 0};

          g_vec3<float_type> pos_0 = {0, 0, 0};
          g_vec3<float_type> pos_1 = {0, 0, 0};
          g_vec3<float_type> pos_2 = {0, 0, 0};
          g_vec3<float_type> pos_3 = {0, 0, 0};

          dvec2 tex_0 = {0.5 * sp_x[i - 1] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j - 1] / sp_y[sp_y.size() - 1]};
          dvec2 tex_1 = {0.5 * sp_x[i - 1] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j] / sp_y[sp_y.size() - 1]};
          dvec2 tex_2 = {0.5 * sp_x[i] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j - 1] / sp_y[sp_y.size() - 1]};
          dvec2 tex_3 = {0.5 * sp_x[i] / sp_x[sp_x.size() - 1], tex_v_shift + 0.35 * sp_y[j] / sp_y[sp_y.size() - 1]};

          g_vec3<float_type> norm_z = {0, 0, 0};
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

          depth[z] = window_depth;
          if (k)
          {
            depth[z] = -window_depth;
          }

          if (i % 2 == 0 && j % 2 == 0 && i != sp_x.size() - 1 && j != sp_y.size() - 1)
          {
            float_type u = 0.2 * 0.5 * (sp_x[i] - sp_x[i - 1]) / sp_x[sp_x.size() - 1];
            float_type v = 0.2 * 0.5 * (sp_y[j] - sp_y[j - 1]) / sp_y[sp_y.size() - 1];
            
            if (!low_quality)
            {
              //alcove
              add_vertex(model, num++, pos_0, norm_x, tex_0 + dvec2{0, v}, only_pos); 
              add_vertex(model, num++, pos_1, norm_x, tex_1 - dvec2{0, v}, only_pos); 
              add_vertex(model, num++, pos_0 + depth, norm_x, tex_0 + dvec2{u, v}, only_pos); 
              add_vertex(model, num++, pos_1 + depth, norm_x, tex_1 + dvec2{u, -v}, only_pos); 
              add_vertex(model, num++, pos_0 + depth, norm_x, tex_0 + dvec2{u, v}, only_pos); 
              add_vertex(model, num++, pos_1, norm_x, tex_1 - dvec2{0, v}, only_pos); 

              add_vertex(model, num++, pos_2, g_vec3<float_type>{0, 0, 0} - norm_x, tex_2 + dvec2{0, v}, only_pos); 
              add_vertex(model, num++, pos_3, g_vec3<float_type>{0, 0, 0} - norm_x, tex_3 - dvec2{0, v}, only_pos); 
              add_vertex(model, num++, pos_2 + depth, g_vec3<float_type>{0, 0, 0} - norm_x, tex_2 + dvec2{-u, v}, only_pos); 
              add_vertex(model, num++, pos_3 + depth, g_vec3<float_type>{0, 0, 0} - norm_x, tex_3 + dvec2{-u, -v}, only_pos); 
              add_vertex(model, num++, pos_2 + depth, g_vec3<float_type>{0, 0, 0} - norm_x, tex_2 + dvec2{-u, v}, only_pos); 
              add_vertex(model, num++, pos_3, g_vec3<float_type>{0, 0, 0} - norm_x, tex_3 - dvec2{0, v}, only_pos);

              add_vertex(model, num++, pos_0, norm_y, tex_0 + dvec2{u, 0}, only_pos); 
              add_vertex(model, num++, pos_2, norm_y, tex_2 - dvec2{u, 0}, only_pos); 
              add_vertex(model, num++, pos_0 + depth, norm_y, tex_0 + dvec2{u, v}, only_pos); 
              add_vertex(model, num++, pos_2 + depth, norm_y, tex_2 + dvec2{-u, v}, only_pos); 
              add_vertex(model, num++, pos_0 + depth, norm_y, tex_0 + dvec2{u, v}, only_pos); 
              add_vertex(model, num++, pos_2, norm_y, tex_2 - dvec2{u, 0}, only_pos); 

              add_vertex(model, num++, pos_1, g_vec3<float_type>{0, 0, 0} - norm_y, tex_1 + dvec2{u, 0}, only_pos); 
              add_vertex(model, num++, pos_3, g_vec3<float_type>{0, 0, 0} - norm_y, tex_3 - dvec2{u, 0}, only_pos); 
              add_vertex(model, num++, pos_1 + depth, g_vec3<float_type>{0, 0, 0} - norm_y, tex_1 + dvec2{u, -v}, only_pos); 
              add_vertex(model, num++, pos_3 + depth, g_vec3<float_type>{0, 0, 0} - norm_y, tex_3 + dvec2{-u, -v}, only_pos); 
              add_vertex(model, num++, pos_1 + depth, g_vec3<float_type>{0, 0, 0} - norm_y, tex_1 + dvec2{u, -v}, only_pos); 
              add_vertex(model, num++, pos_3, g_vec3<float_type>{0, 0, 0} - norm_y, tex_3 - dvec2{u, 0}, only_pos);



              //window
              int wn = windows.size()/FLOAT_PER_VERTEX;
              add_vertex(windows, wn++, pos_0 + depth, norm_z, tex_0, only_pos); 
              add_vertex(windows, wn++, pos_1 + depth, norm_z, tex_1, only_pos); 
              add_vertex(windows, wn++, pos_2 + depth, norm_z, tex_2, only_pos); 
              add_vertex(windows, wn++, pos_3 + depth, norm_z, tex_3, only_pos); 
              add_vertex(windows, wn++, pos_2 + depth, norm_z, tex_2, only_pos); 
              add_vertex(windows, wn++, pos_1 + depth, norm_z, tex_1, only_pos); 
            }
            tex_0 += dvec2{u, v};
            tex_1 += dvec2{u, -v};
            tex_2 += dvec2{-u, v};
            tex_3 += dvec2{-u, -v};
            
          }
          else
          {
            //facade segment
            add_vertex(model, num++, pos_0, norm_z, tex_0, only_pos); 
            add_vertex(model, num++, pos_1, norm_z, tex_1, only_pos); 
            add_vertex(model, num++, pos_2, norm_z, tex_2, only_pos); 
            add_vertex(model, num++, pos_3, norm_z, tex_3, only_pos); 
            add_vertex(model, num++, pos_2, norm_z, tex_2, only_pos); 
            add_vertex(model, num++, pos_1, norm_z, tex_1, only_pos); 

            //inner side of the wall
            if (!low_quality)
            {
              add_vertex(model, num++, pos_0 + depth, -norm_z, tex_0, only_pos); 
              add_vertex(model, num++, pos_1 + depth, -norm_z, tex_1, only_pos); 
              add_vertex(model, num++, pos_2 + depth, -norm_z, tex_2, only_pos); 
              add_vertex(model, num++, pos_3 + depth, -norm_z, tex_3, only_pos); 
              add_vertex(model, num++, pos_2 + depth, -norm_z, tex_2, only_pos); 
              add_vertex(model, num++, pos_1 + depth, -norm_z, tex_1, only_pos); 
            }
          }
        }
      }
    }
  }
  template void get_walls_from_splines<dfloat>(std::vector<dfloat> &model, std::vector<dfloat> &windows,
                              const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, 
                              dfloat length_z, dfloat window_depth, int x, int y, int z, float tex_v_shift, 
                              bool low_quality, bool only_pos);

  template<typename float_type>
  inline void create_floor_simple(std::vector<float_type> &model, const std::vector<float_type> &sp_x, const std::vector<float_type> &sp_y, 
                                  const std::vector<float_type> &sp_z, bool only_pos)
  {
    int num = model.size()/FLOAT_PER_VERTEX;
    add_vertex(model, num++, g_vec3<float_type>{0, 0, 0}, g_vec3<float_type>{0, -1, 0}, dvec2{0, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{0, 0, sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, -1, 0}, dvec2{0, 1}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], 0, 0}, g_vec3<float_type>{0, -1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], 0, sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, -1, 0}, dvec2{0.5, 1}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], 0, 0}, g_vec3<float_type>{0, -1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{0, 0, sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, -1, 0}, dvec2{0, 1}, only_pos); 
  }
  template void create_floor_simple<dfloat>(std::vector<dfloat> &model, const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, 
                                  const std::vector<dfloat> &sp_z, bool only_pos);

  template<typename float_type>
  inline void create_roof_simple(std::vector<float_type> &model, const std::vector<float_type> &sp_x, const std::vector<float_type> &sp_y, 
                                 const std::vector<float_type> &sp_z, bool only_pos)
  {
    int num = model.size()/FLOAT_PER_VERTEX;
    add_vertex(model, num++, g_vec3<float_type>{0, sp_y[sp_y.size() - 1], 0}, g_vec3<float_type>{0, 1, 0}, dvec2{0.5, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{0, sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, 1, 0}, dvec2{0.5, 1}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], 0}, g_vec3<float_type>{0, 1, 0}, dvec2{1, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, 1, 0}, dvec2{1, 1}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{sp_x[sp_x.size() - 1], sp_y[sp_y.size() - 1], 0}, g_vec3<float_type>{0, 1, 0}, dvec2{1, 0.7}, only_pos); 
    add_vertex(model, num++, g_vec3<float_type>{0, sp_y[sp_y.size() - 1], sp_z[sp_z.size() - 1]}, g_vec3<float_type>{0, 1, 0}, dvec2{0.5, 1}, only_pos);
  }
  template void create_roof_simple<dfloat>(std::vector<dfloat> &model, const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, 
                                 const std::vector<dfloat> &sp_z, bool only_pos);

  template<typename float_type>
  void splines_to_building(std::vector<float_type> &model, std::vector<float_type> &windows,
                           const std::vector<float_type> &sp_x, const std::vector<float_type> &sp_y, 
                           const std::vector<float_type> &sp_z, float_type window_depth, bool low_quality, bool only_pos)
  {
    get_walls_from_splines(model, windows, sp_x, sp_y, sp_z[sp_z.size() - 1], window_depth, 0, 1, 2, 0, low_quality, only_pos);
    get_walls_from_splines(model, windows, sp_z, sp_y, sp_x[sp_x.size() - 1], window_depth, 2, 1, 0, 0.35, low_quality, only_pos);
    
    create_floor_simple(model, sp_x, sp_y, sp_z, only_pos);
    create_roof_simple(model, sp_x, sp_y, sp_z, only_pos);
  }
  template void splines_to_building<dfloat>(std::vector<dfloat> &model, std::vector<dfloat> &windows,
                           const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, 
                           const std::vector<dfloat> &sp_z, dfloat window_depth, bool low_quality, bool only_pos);

  template<typename float_type>
  void splines_to_box(std::vector<float_type> &model, const std::vector<float_type> &sp_x, const std::vector<float_type> &sp_y, 
                           const std::vector<float_type> &sp_z, bool only_pos)
  {
    auto add_quad = [&only_pos, &model](const g_vec3<float_type> &p, const g_vec3<float_type> &v1, const g_vec3<float_type> &v2, const g_vec3<float_type> &n)
    {
      int num = model.size()/FLOAT_PER_VERTEX;
      add_vertex(model, num++, p          , n, dvec2(0,0), only_pos); 
      add_vertex(model, num++, p + v1     , n, dvec2(0,0), only_pos); 
      add_vertex(model, num++, p + v2     , n, dvec2(0,0), only_pos); 
      add_vertex(model, num++, p + v1 + v2, n, dvec2(0,0), only_pos); 
      add_vertex(model, num++, p + v2     , n, dvec2(0,0), only_pos); 
      add_vertex(model, num++, p + v1     , n, dvec2(0,0), only_pos); 
    };

    add_quad(g_vec3<float_type>(0,0,0), g_vec3<float_type>(sp_x.back(), 0, 0), g_vec3<float_type>(0, sp_y.back(), 0), g_vec3<float_type>(0,0,-1));
    add_quad(g_vec3<float_type>(0,0,sp_z.back()), g_vec3<float_type>(0, sp_y.back(), 0), g_vec3<float_type>(sp_x.back(), 0, 0), g_vec3<float_type>(0,0,1));

    add_quad(g_vec3<float_type>(0,0,0), g_vec3<float_type>(0, sp_y.back(), 0), g_vec3<float_type>(0, 0, sp_z.back()), g_vec3<float_type>(-1,0,0));
    add_quad(g_vec3<float_type>(sp_x.back(),0,0), g_vec3<float_type>(0, 0, sp_z.back()), g_vec3<float_type>(0, sp_y.back(), 0), g_vec3<float_type>(1,0,0));
    create_floor_simple(model, sp_x, sp_y, sp_z, only_pos);
    create_roof_simple(model, sp_x, sp_y, sp_z, only_pos);
  }
  template void splines_to_box<dfloat>(std::vector<dfloat> &model, const std::vector<dfloat> &sp_x, const std::vector<dfloat> &sp_y, 
                           const std::vector<dfloat> &sp_z, bool only_pos);

  template<typename float_type>
  PartOffsets create_building(const std::vector<float_type> &params, std::vector<float_type> &vert, ModelQuality quality)
  {
    enum BuildingQuality
    {
      BQ_BOXES,
      BQ_FACADES,
      BQ_FULL,
      BQ_COUNT
    };
    BuildingQuality bq = (BuildingQuality)CLAMP(quality.quality_level, 0, BQ_COUNT-1);
    std::vector<float_type> spline_x = create_wall_spline(params[0], params[1], params[2], params[3], params[4]);
    std::vector<float_type> spline_y = create_wall_spline(params[5], params[6], params[7], params[8], params[9]);
    std::vector<float_type> spline_z = create_wall_spline(params[10], params[11], params[12], params[13], params[14]);

    std::vector<float_type> windows_mesh;
    
    if (bq == BQ_BOXES)
      splines_to_box(vert, spline_x, spline_y, spline_z, quality.create_only_position);
    else
      splines_to_building(vert, windows_mesh, spline_x, spline_y, spline_z, params[15], bq == BQ_FACADES, quality.create_only_position);

    int windows_offset = vert.size();
    if (bq != BQ_BOXES)
    {
      vert.reserve(vert.size() + windows_mesh.size());
      for (auto &v : windows_mesh)
        vert.push_back(v);
    }

    float_type total_size_x = params[0] + params[4]*(params[1] + params[3]) + params[2];
    float_type total_size_y = params[5] + params[9]*(params[6] + params[8]) + params[7];
    float_type total_size_z = params[10] + params[14]*(params[11] + params[13]) + params[12];
    g_mat43<float_type> sc2 = scale(translate(ident<float_type>(), g_vec3<float_type>{-0.5, 0, 0}), g_vec3<float_type>{1.0/total_size_x, params[16]/total_size_y, params[17]/total_size_z});
    transform(vert, sc2);

    return {
            {"main_part", 0},
            {"windows", windows_offset}
           };
  }
  template PartOffsets create_building<dfloat>(const std::vector<dfloat> &params, std::vector<dfloat> &vert, ModelQuality quality);
};