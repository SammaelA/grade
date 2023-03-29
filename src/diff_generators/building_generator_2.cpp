#include "building_generator_2.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  #define vec2 g_vec2<real>
  #define vec3 g_vec3<real>
  #define vec4 g_vec4<real>
  #define mat43 g_mat43<real>

  template <typename real>
  PartOffsets create_building_2t(const std::vector<real> &params, std::vector<real> &model, ModelQuality quality)
  {
    struct Quad
    {
      Quad() = default;
      Quad(const Quad &) = default;
      Quad(const vec3 &_p1, const vec2 &_tc_1,
           const vec3 &_v1, const vec2 &_tc_v1,
           const vec3 &_v2, const vec2 &_tc_v2,
           const vec3 &_n) : p1(_p1), v1(_v1), v2(_v2), tc(_tc_1), tc_v1(_tc_v1), tc_v2(_tc_v2)
      {
        n = normalize(_n);
      }
      Quad(const vec3 &_p1,
           const vec3 &_v1,
           const vec3 &_v2,
           const vec3 &_n) : p1(_p1), v1(_v1), v2(_v2), tc(vec2(0,0)), tc_v1(vec2(1,0)), tc_v2(vec2(0,1))
      {
        n = normalize(_n);
      }
    Quad &operator=(const Quad &) = default;
    Quad &operator=(Quad &&) = default;
      vec3 p1, v1, v2, n;
      vec2 tc, tc_v1, tc_v2;
    };

    auto add_vertex = [&](std::vector<real> &M, const vec3 &pos, const vec3 &n, const vec2 &tc)
    {
      //we are rely heavily on reserve() on vectors and assume, that most push_back()
      //functions will not lead to increasing vector size and reallocation

      M.push_back(pos.x);
      M.push_back(pos.y);
      M.push_back(pos.z);

      M.push_back(n.x);
      M.push_back(n.y);
      M.push_back(n.z);

      M.push_back(tc.x);
      M.push_back(tc.y);
    };

    auto make_quad = [&](std::vector<real> &M, const Quad &q)
    {
      /*
      p1+v2***p1+v1+v2
      *       *
      *       *
      *       *
      *       *
      *       *
      p1 *****p1+v1
      */

      add_vertex(M, q.p1, q.n, q.tc);
      add_vertex(M, q.p1+q.v1, q.n, q.tc+q.tc_v1);
      add_vertex(M, q.p1+q.v2, q.n, q.tc+q.tc_v2);
      add_vertex(M, q.p1+q.v1, q.n, q.tc+q.tc_v1);
      add_vertex(M, q.p1+q.v1+q.v2, q.n, q.tc+q.tc_v1+q.tc_v2);
      add_vertex(M, q.p1+q.v2, q.n, q.tc+q.tc_v2);
    };

    auto make_panel = [&](std::vector<real> &M_ext, std::vector<real> &M_int, 
                          bool has_interior, real thick,
                          const Quad &q)
    {
      make_quad(M_ext, q);
      if (has_interior)
      {
        Quad q2(q.p1 - thick*q.n, q.tc, q.v2, q.tc_v2, q.v1, q.tc_v1, -q.n);
        make_quad(M_int, q2);
      }
    };

    auto make_box = [&](std::vector<real> &M, real thick, const Quad &q)
    {
      //Quad MUST be a rectangle (dot(v1, v2) = dot(tc_v1, tc_v2) = 0)
      //tc mentioned in quad will be used to cover ALL the box on texture
      //without stretching (not all texture will be used consequently)
      //better to have ||v1|| >= ||v2|| >= thick, but not mandatory

      real a = length(q.v1);
      real b = length(q.v2);
      real c = thick;

      vec2 tc_v1 = q.tc_v1 / (2*(b+c));
      vec2 tc_v2 = q.tc_v2 / (a + 2*c);

      Quad side_l = Quad(q.p1, q.tc + c*tc_v1 + c*tc_v2, q.v2, a*tc_v2, -thick*q.n, -c*tc_v1, -normalize(q.v1));
      Quad front = Quad(q.p1, q.tc + c*tc_v1 + c*tc_v2, q.v1, b*tc_v1, q.v2, a*tc_v2, q.n);
      Quad side_r = Quad(q.p1 + q.v1, q.tc + (b+c)*tc_v1 + c*tc_v2, -thick*q.n, c*tc_v1, q.v2, a*tc_v2, normalize(q.v1));
      Quad back = Quad(q.p1 - thick*q.n, q.tc + (b+2*c)*tc_v1 + c*tc_v2, q.v2, a*tc_v2, q.v1, b*tc_v1, -q.n);
      Quad bottom = Quad(q.p1, q.tc + c*tc_v1 + c*tc_v2, -thick*q.n, -c*tc_v2, q.v1, b*tc_v1, -normalize(q.v2));
      Quad top = Quad(q.p1 + q.v2, q.tc + c*tc_v1 + (a+c)*tc_v2, q.v1, b*tc_v1, -thick*q.n, c*tc_v2, normalize(q.v2));

      make_quad(M, side_l);
      make_quad(M, front);
      make_quad(M, side_r);
      make_quad(M, back);
      make_quad(M, bottom);
      make_quad(M, top);
    };

    enum WindowQuality
    {
      WQ_LOW, //only glass quad
      WQ_HIGH //glass quad + alcove + frame
    };

    enum WindowSplit
    {
      WS_LEFT,
      WS_RIGHT,
      WS_NONE,
      WS_BOTH
    };

    enum BalconyQuality
    {
      BQ_NONE,//no balcony at all
      BQ_LOW, //plates instead of bars and windows
      BQ_HIGH //full quality
    };

    auto make_window = [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                           std::vector<real> &M_frame,
                           WindowQuality wq, real hor_parts, real vert_parts, 
                           WindowSplit ws, real split_h_q, real wall_thick, 
                           real outer_frame_width, real inner_frame_width,
                           real glass_deep_q, const Quad &q)
    {
      if (wq == WindowQuality::WQ_LOW)
      {
        make_quad(M_glass, q);
      }
      else
      {
        Quad deepened_q = q;
        deepened_q.p1 -= glass_deep_q*wall_thick*q.n;
        make_quad(M_glass, deepened_q);

        // outer frame
        real w = outer_frame_width * length(q.v2) / length(q.v1);
        real h = outer_frame_width;
        real th = 2*glass_deep_q*wall_thick;
        vec3 p_off = q.p1;

        make_box(M_frame, th, Quad(p_off, vec2(0, 0), w * q.v1, vec2(1, 0), q.v2, vec2(0, 1), q.n));
        make_box(M_frame, th, Quad(p_off + (1 - w) * q.v1, vec2(0, 0), w * q.v1, vec2(1, 0), q.v2, vec2(0, 1), q.n));
        make_box(M_frame, th, Quad(p_off + w * q.v1, vec2(0, 0), (1 - 2 * w) * q.v1, vec2(1, 0), h * q.v2, vec2(0, 1), q.n));
        make_box(M_frame, th, Quad(p_off + w * q.v1 + (1 - h) * q.v2, vec2(0, 0), (1 - 2 * w) * q.v1, vec2(1, 0), h * q.v2, vec2(0, 1), q.n));

        // inner frames
        real w1 = inner_frame_width * length(q.v2) / length(q.v1);
        real h1 = inner_frame_width;
        for (int i = 1; i < hor_parts; i++)
          make_box(M_frame, th, Quad(p_off + w * q.v1 + (i/hor_parts - 0.5*h1) * q.v2, vec2(0, 0), (1 - 2 * w) * q.v1, vec2(1, 0), h1 * q.v2, vec2(0, 1), q.n));
        for (int i = 1; i < vert_parts; i++)
          make_box(M_frame, th, Quad(p_off + h * q.v2 + (i/vert_parts - 0.5*w1) * q.v1, vec2(0, 0), w1 * q.v1, vec2(1, 0), (1 - 2 * h) * q.v2, vec2(0, 1), q.n));

        if (ws == WindowSplit::WS_LEFT || ws == WindowSplit::WS_BOTH)
        {
          make_box(M_frame, th, Quad(p_off + w * q.v1 + split_h_q * q.v2, vec2(0, 0), (1/vert_parts - 0.5*w1 - w) * q.v1, vec2(1, 0), h1 * q.v2, vec2(0, 1), q.n));
        }

        if (ws == WindowSplit::WS_RIGHT || ws == WindowSplit::WS_BOTH)
        {
          make_box(M_frame, th, Quad(p_off + (1 - 1/vert_parts + 0.5*w1) * q.v1 + split_h_q * q.v2, vec2(0, 0), (1/vert_parts - 0.5*w1 - w) * q.v1, vec2(1, 0), h1 * q.v2, vec2(0, 1), q.n));
        }

        //alcove
        make_quad(M_wall, Quad(q.p1, vec2(0, 0), q.v1, vec2(1, 0), -wall_thick*q.n, vec2(0, 1), q.v2));
        make_quad(M_wall, Quad(q.p1 + q.v2, vec2(0, 0), -wall_thick*q.n, vec2(0, 1), q.v1, vec2(1, 0), -q.v2));
        make_quad(M_wall, Quad(q.p1, vec2(0, 0), -wall_thick*q.n, vec2(0, 1), q.v2, vec2(1, 0), q.v1));
        make_quad(M_wall, Quad(q.p1 + q.v1, vec2(0, 0), q.v2, vec2(1, 0), -wall_thick*q.n, vec2(0, 1), -q.v1));
      }
    };

    auto make_stripe_wall = [&](std::vector<real> &M_ext, std::vector<real> &M_int, 
                                bool has_interior, real thick, real bottom_offset_q,
                                real top_offset_q, real floors_count,
                                const Quad &q)
    {
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      Quad q_part = q;

      //bottom part
      real t = 0;
      for (int i=0;i<bottom_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_ext, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (bottom_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (bottom_offset_q - t)*q.tc_v2;
      make_panel(M_ext, M_int, has_interior, thick, q_part);

      //floors
      q_part.p1 = q.p1 + bottom_offset_q * y_mul * q.v2;
      q_part.v2 = y_mul * q.v2;
      q_part.tc_v2 = q.tc_v2;

      for (int i=0;i<floors_count;i++)
      {
        make_panel(M_ext, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
      }

      //top part
      t = 0;
      for (int i=0;i<top_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_ext, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (top_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (top_offset_q - t)*q.tc_v2;
      make_panel(M_ext, M_int, has_interior, thick, q_part);
    };

    auto make_window_panel = [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                                 std::vector<real> &M_frame, std::vector<real> &M_int,
                                 WindowQuality wq, real hor_parts, real vert_parts,
                                 WindowSplit ws, real split_h_q, real wall_thick,
                                 real outer_frame_width, real inner_frame_width,
                                 real glass_deep_q,
                                 real window_width_q, real window_height_q,
                                 const Quad &q)
    {
      real w = 0.5f * (1 - window_width_q) * length(q.v2) / length(q.v1);
      real h = 0.5f * (1 - window_height_q);

      bool int_wall = (wq == WindowQuality::WQ_HIGH);

      make_panel(M_wall, M_int, int_wall, wall_thick, Quad(q.p1, vec2(0, 0), w * q.v1, vec2(w, 0), q.v2, vec2(0, 1), q.n));
      make_panel(M_wall, M_int, int_wall, wall_thick, Quad(q.p1 + (1 - w) * q.v1, vec2(1-w, 0), w * q.v1, vec2(w, 0), q.v2, vec2(0, 1), q.n));
      make_panel(M_wall, M_int, int_wall, wall_thick, Quad(q.p1 + w * q.v1, vec2(w, 0), (1 - 2 * w) * q.v1, vec2((1 - 2 * w), 0), h * q.v2, vec2(0, h), q.n));
      make_panel(M_wall, M_int, int_wall, wall_thick, Quad(q.p1 + w * q.v1 + (1 - h) * q.v2, vec2(w, 1-h), (1 - 2 * w) * q.v1, vec2((1 - 2 * w), 0), h * q.v2, vec2(0, h), q.n));

      Quad window_q(q.p1 + w * q.v1 + h * q.v2, vec2(0, 0), (1 - 2 * w) * q.v1, vec2(1, 0), (1 - 2 * h) * q.v2, vec2(0, 1), q.n);
      make_window(M_glass, M_int, M_frame, wq, hor_parts, vert_parts, ws, split_h_q, wall_thick, outer_frame_width, inner_frame_width, glass_deep_q, window_q);
    };

    auto make_stripe_windows = [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                                 std::vector<real> &M_frame, std::vector<real> &M_int,
                                 WindowQuality wq, real hor_parts, real vert_parts,
                                 WindowSplit ws, real split_h_q, 
                                 real outer_frame_width, real inner_frame_width,
                                 real glass_deep_q,
                                 real window_width_q, real window_height_q,
                                 real thick, real bottom_offset_q, real top_offset_q, real floors_count,
                                 const Quad &q)
    {
      bool has_interior = (wq == WindowQuality::WQ_HIGH);
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      Quad q_part = q;

      //bottom part
      real t = 0;
      for (int i=0;i<bottom_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_wall, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (bottom_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (bottom_offset_q - t)*q.tc_v2;
      make_panel(M_wall, M_int, has_interior, thick, q_part);

      //floors
      q_part.p1 = q.p1 + bottom_offset_q * y_mul * q.v2;
      q_part.v2 = y_mul * q.v2;
      q_part.tc_v2 = q.tc_v2;

      for (int i=0;i<floors_count;i++)
      {
        make_window_panel(M_glass, M_wall, M_frame, M_int, wq, hor_parts, vert_parts, ws, split_h_q, thick, outer_frame_width, inner_frame_width, glass_deep_q, 
                          window_width_q, window_height_q, q_part);
        q_part.p1 += y_mul*q.v2;
      }

      //top part
      t = 0;
      for (int i=0;i<top_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_wall, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (top_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (top_offset_q - t)*q.tc_v2;
      make_panel(M_wall, M_int, has_interior, thick, q_part);
    };

    auto make_door = [&](std::vector<real> &M_door, real thick, const Quad &q)
    {
      make_box(M_door, thick, q);
    };

    auto make_stripe_entrance = [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                                    std::vector<real> &M_frame, std::vector<real> &M_int,
                                    std::vector<real> &M_door,
                                    WindowQuality wq, real hor_parts, real vert_parts,
                                    WindowSplit ws, real split_h_q,
                                    real outer_frame_width, real inner_frame_width,
                                    real glass_deep_q,
                                    real window_width_q, real window_height_q,
                                    real thick, real bottom_offset_q, real top_offset_q, real floors_count,
                                    real door_offset_q, real door_width_q, real door_height_q, real door_deep_q,
                                    real stairs_h_q, real stairs_w_q,
                                    const Quad &q)
    {
      bool has_interior = (wq == WindowQuality::WQ_HIGH);
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      Quad q_part = q;

      //door
      
      real w_off = 0.5f * (1 - door_width_q * length(y_mul*q.v2) / length(q.v1));
      real door_wq = door_width_q * length(y_mul*q.v2) / length(q.v1);
      Quad q_door(q.p1 + door_offset_q*y_mul*q.v2 + w_off*q.v1 - door_deep_q*thick*q.n, door_wq * q.v1, door_height_q*y_mul*q.v2,q.n);
      make_door(M_door, 0.2*thick, q_door);
      
      //doorway
      {
        Quad q_l(q.p1 + w_off*q.v1, -thick*q.n, (bottom_offset_q+1) * y_mul * q.v2, q.v1);
        Quad q_r(q.p1 + (w_off + door_wq)*q.v1, (bottom_offset_q+1) * y_mul * q.v2, -thick*q.n, -q.v1);
        make_quad(M_int, q_l);
        make_quad(M_int, q_r);
      }

      //walls left and right from the door
      real t = 0;
      for (int i=0;i<bottom_offset_q+1 - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        {
          Quad q_part1 = q_part;
          q_part1.v1 = w_off*q.v1;
          q_part1.tc_v1 = vec2(w_off,0);

          Quad q_part2 = q_part;
          q_part2.v1 = w_off*q.v1;
          q_part2.p1 += (w_off + door_wq)*q.v1;
          q_part2.tc_v1 = vec2(w_off,0);
          q_part2.tc = vec2(1 - w_off,0);

          make_panel(M_wall, M_int, has_interior, thick, q_part1);
          make_panel(M_wall, M_int, has_interior, thick, q_part2);
        }
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (bottom_offset_q+1 - t)*y_mul*q.v2;
      q_part.tc_v2 = (bottom_offset_q+1 - t)*q.tc_v2;
      {
        Quad q_part1 = q_part;
        q_part1.v1 = w_off * q.v1;
        q_part1.tc_v1 = vec2(w_off, 0);
        Quad q_part2 = q_part;
        q_part2.v1 = w_off * q.v1;
        q_part2.p1 += (w_off + door_wq) * q.v1;
        q_part2.tc_v1 = vec2(w_off, 0);
        q_part2.tc = vec2(1 - w_off, 0);
        make_panel(M_wall, M_int, has_interior, thick, q_part1);
        make_panel(M_wall, M_int, has_interior, thick, q_part2);
      }

      //wall above the door
      {
        Quad q_above(q.p1 + (door_offset_q+door_height_q)*y_mul*q.v2 + w_off*q.v1, vec2(w_off, 0), door_wq * q.v1, vec2(1-2*w_off,0), 
                     (bottom_offset_q+1 - door_height_q - door_offset_q)*y_mul*q.v2, vec2(0,bottom_offset_q+1 - door_height_q),q.n);
        make_panel(M_wall, M_int, has_interior, thick, q_above);
      
        Quad q_doorway_above(q.p1 + (door_offset_q+door_height_q)*y_mul*q.v2 + w_off*q.v1, -thick*q.n, door_wq * q.v1, -q.v2);
        make_quad(M_int, q_doorway_above);
      }

      //box under the door and stairs
      {
        Quad q_under(q.p1 + w_off*q.v1, door_wq * q.v1, door_offset_q * y_mul * q.v2, q.n);
        make_box(M_int, thick, q_under);

        real sw = stairs_w_q * length(y_mul*q.v2);
        Quad s_quad(q.p1 + w_off*q.v1, door_wq * q.v1, q.v2, q.n);
        real sh = door_offset_q;
        while (sh > 0)
        {
          s_quad.p1 += sw*q.n;
          s_quad.v2 = sh*y_mul*q.v2;
          make_box(M_int, sw, s_quad);

          sh -= stairs_h_q;
        }
      }

      //floors
      q_part.p1 = q.p1 + (bottom_offset_q+1) * y_mul * q.v2;
      q_part.v2 = y_mul * q.v2;
      q_part.tc_v2 = q.tc_v2;

      for (int i=1;i<floors_count;i++)
      {
        make_window_panel(M_glass, M_wall, M_frame, M_int, wq, hor_parts, vert_parts, ws, split_h_q, thick, outer_frame_width, inner_frame_width, glass_deep_q, 
                            window_width_q, window_height_q, q_part);
        q_part.p1 += y_mul*q.v2;
      }

      //top part
      t = 0;
      for (int i=0;i<top_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_wall, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (top_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (top_offset_q - t)*q.tc_v2;
      make_panel(M_wall, M_int, has_interior, thick, q_part);
    };

    auto make_balcony = [&](std::vector<real> &M_plate, std::vector<real> &M_bars, BalconyQuality bq,
                            real plate_thick, real depth, real bars_thick, real bars_distance, real bars_height,
                            const Quad &q)
    {
      if (bq > BalconyQuality::BQ_NONE)
      {
        Quad q_plate(q.p1 + depth*q.n, q.v1, -depth*q.n, q.v2);
        make_box(M_plate, plate_thick, q_plate);

        if (bq > BalconyQuality::BQ_LOW)
        {
          real t = 0;
          real b_th = bars_thick*length(q.v1);
          real b_th_y = bars_thick*length(q.v1) / length(q.v2);
          //bars in front
          while (t < 1)
          {
            make_box(M_bars, b_th, Quad(q.p1 + depth*q.n + t*q.v1, bars_thick*q.v1, bars_height*q.v2, q.n));
            t += bars_distance;
          }
          make_box(M_bars, b_th, Quad(q.p1 + depth*q.n + bars_height*q.v2, q.v1, b_th_y*q.v2, q.n));

          t = 0;
          real b_dist = bars_distance*length(q.v1);
          while (t < depth)
          {
            make_box(M_bars, b_th, Quad(q.p1 + (depth-t)*q.n, bars_thick*q.v1, bars_height*q.v2, q.n));
            make_box(M_bars, b_th, Quad(q.p1 + (depth-t)*q.n + (1-bars_thick)*q.v1, bars_thick*q.v1, bars_height*q.v2, q.n));
            t+= b_dist;
          }
          make_box(M_bars, depth - b_th, Quad(q.p1 + (depth-b_th)*q.n + bars_height*q.v2, bars_thick*q.v1, b_th_y*q.v2, q.n));
          make_box(M_bars, depth - b_th, Quad(q.p1 + (depth-b_th)*q.n + (1-bars_thick)*q.v1 + bars_height*q.v2, bars_thick*q.v1, b_th_y*q.v2, q.n));
        }
        else
        {

        }
      }
    };

    auto make_stripe_balcony =  [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                                    std::vector<real> &M_frame, std::vector<real> &M_int,
                                    std::vector<real> &M_bars,
                                    WindowQuality wq, real hor_parts, real vert_parts,
                                    WindowSplit ws, real split_h_q,
                                    real outer_frame_width, real inner_frame_width,
                                    real glass_deep_q,
                                    real window_width_q, real window_height_q,
                                    real thick, real bottom_offset_q, real top_offset_q,
                                    BalconyQuality bq,
                                    real floors_count, real balcony_height_q, real balcony_depth_q,
                                    real bars_thick, real bars_distance, real balcony_start_floor,
                                    const Quad &q)
    {
      bool has_interior = (wq == WindowQuality::WQ_HIGH);
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      Quad q_part = q;

      //bottom part
      real t = 0;
      for (int i=0;i<bottom_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_wall, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (bottom_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (bottom_offset_q - t)*q.tc_v2;
      make_panel(M_wall, M_int, has_interior, thick, q_part);

      //floors
      q_part.p1 = q.p1 + bottom_offset_q * y_mul * q.v2;
      q_part.v2 = y_mul * q.v2;
      q_part.tc_v2 = q.tc_v2;

      for (int i=0;i<floors_count;i++)
      {
        //make window
        make_window_panel(M_glass, M_wall, M_frame, M_int, wq, hor_parts, vert_parts, ws, split_h_q, thick, outer_frame_width, inner_frame_width, glass_deep_q, 
                          window_width_q, window_height_q, q_part);


        //make balcony
        if (i >= balcony_start_floor)
          make_balcony(M_int, M_bars, bq, 0.04 * length(y_mul*q.v2), balcony_depth_q * length(y_mul*q.v2), bars_thick, bars_distance, balcony_height_q, q_part);

        q_part.p1 += y_mul*q.v2;
      }

      //top part
      t = 0;
      for (int i=0;i<top_offset_q - 1;i++)
      {
        q_part.v2 = y_mul * q.v2;
        make_panel(M_wall, M_int, has_interior, thick, q_part);
        q_part.p1 += y_mul*q.v2;
        t = i + 1;
      }
      q_part.v2 = (top_offset_q - t)*y_mul*q.v2;
      q_part.tc_v2 = (top_offset_q - t)*q.tc_v2;
      make_panel(M_wall, M_int, has_interior, thick, q_part);
    };

    auto make_roof_segment = [&](std::vector<real> &M_wall, std::vector<real> &M_roof, 
                                 real roof_base_h, real roof_angle_rad, 
                                 bool left_side, bool right_side,
                                 const Quad &q)
    {
      Quad q_rp = q;
      q_rp.p1 += roof_base_h*q.n;
      make_box(M_wall, length(q.v2), Quad(q.p1, q.v1, roof_base_h*q.n, -q.v2));

      real l1 = 0.5f*length(q.v2);
      real h1 = l1*tan(roof_angle_rad);
      vec3 v1 = ((real)0.5f)*q.v2 + h1*q.n;
      make_quad(M_roof, Quad(q.p1 + roof_base_h*q.n, q.v1, v1, cross(q.v1, v1)));
      v1 = ((real)-0.5f)*q.v2 + h1*q.n;
      make_quad(M_roof, Quad(q.p1 + roof_base_h*q.n + q.v2, q.v1, v1, -cross(q.v1, v1)));

      if (left_side)
      {
        add_vertex(M_wall, q.p1 + roof_base_h*q.n, -q.v1, vec2(0,0));
        add_vertex(M_wall, q.p1 + (roof_base_h)*q.n + q.v2 + v1, -q.v1, vec2(0, tan(roof_angle_rad)));
        add_vertex(M_wall, q.p1 + roof_base_h*q.n + q.v2, -q.v1, vec2(1,0));
      }

      if (right_side)
      {
        add_vertex(M_wall, q.p1 + q.v1 + roof_base_h*q.n, q.v1, vec2(0,0));
        add_vertex(M_wall, q.p1 + q.v1 + (roof_base_h)*q.n + q.v2 + v1, q.v1, vec2(0, tan(roof_angle_rad)));
        add_vertex(M_wall, q.p1 + q.v1 + roof_base_h*q.n + q.v2, q.v1, vec2(1,0));
      }
    };

    auto make_insides = [&](std::vector<real> &M_int, real depth, bool left_wall,
                            real thick, real bottom_offset_q, real top_offset_q, real floors_count,
                            const Quad &q /*stripe quad*/)
    {
      //floors (and roofs)
      vec3 th = 0.1*thick*normalize(q.v2);
      vec3 p1 = q.p1 + 0.1*thick*normalize(q.v1);
      vec3 v1 = q.v1 - 0.1*thick*normalize(q.v1);
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      make_box(M_int, depth - 2*thick, Quad(p1 - q.n*thick, v1, th, q.n));
      p1 += y_mul*bottom_offset_q*q.v2;

      for (int i=0;i<=floors_count;i++)
      {
        make_box(M_int, depth - 2*thick, Quad(p1 - q.n*thick, v1, th, q.n));
        p1 += y_mul*q.v2;
      }

      //middle wall
      make_box(M_int, 0.1*thick, Quad(q.p1 - 0.5*depth*q.n, v1, q.v2, q.n));

      if (left_wall)
      {
        make_box(M_int, depth - 2*thick, Quad(q.p1 - q.n*thick, 0.1*thick*normalize(q.v1), q.v2, q.n));
      }
    };

    enum ParameterNames
    {
      I_ENTRANCES_COUNT,
      I_FLOORS_COUNT,
      F_WALL_THICKNESS,
      F_BOTTOM_OFFSET_Q,
      F_TOP_OFFSET_Q,

      I_ES_STRIPES, //ES referes for End Section
      I_ES_CODE,
      F_ES_WALL_STRIPE_SIZE,

      I_MS_STRIPES, //MS referes for Middle Section
      I_MS_CODE,
      F_MS_WALL_STRIPE_SIZE,

      I_SS_STRIPES, //SS referes for Side Section
      I_SS_CODE,
      F_SS_WALL_STRIPE_SIZE,

      F_WINDOW_OUTER_FRAME_WIDTH,
      F_WINDOW_INNER_FRAME_WIDTH,
      F_WINDOW_GLASS_DEPTH,

      F_ROOF_BASE_H,
      F_ROOF_SLOPE,
      F_ROOF_OVERSIZE,

      I_BASE_WINDOW_HOR_SPLITS,
      I_BASE_WINDOW_VERT_SPLITS,
      I_BASE_WINDOW_SPLIT_TYPE,
      F_BASE_WINDOW_SPLIT_H,
      F_BASE_WINDOW_WIDTH_Q,
      F_BASE_WINDOW_HEIGHT_Q,
      F_BASE_WINDOW_STRIPE_SIZE,

      I_BALCONY_WINDOW_HOR_SPLITS,
      I_BALCONY_WINDOW_VERT_SPLITS,
      I_BALCONY_WINDOW_SPLIT_TYPE,
      F_BALCONY_WINDOW_SPLIT_H,
      F_BALCONY_WINDOW_WIDTH_Q,
      F_BALCONY_WINDOW_HEIGHT_Q,
      F_BALCONY_WINDOW_STRIPE_SIZE,
      F_BALCONY_HEIGHT_Q,
      F_BALCONY_DEPTH_Q,
      F_BALCONY_BARS_THICKNESS,
      F_BALCONY_BARS_DISTANCE,
      I_BALCONY_START_FLOOR,

      I_ENTRANCE_WINDOW_HOR_SPLITS,
      I_ENTRANCE_WINDOW_VERT_SPLITS,
      I_ENTRANCE_WINDOW_SPLIT_TYPE,
      F_ENTRANCE_WINDOW_SPLIT_H,
      F_ENTRANCE_WINDOW_WIDTH_Q,
      F_ENTRANCE_WINDOW_HEIGHT_Q,
      F_ENTRANCE_WINDOW_STRIPE_SIZE,
      F_DOOR_OFFSET_Q,
      F_DOOR_WIDTH_Q,
      F_DOOR_HEIGHT_Q,
      F_DOOR_DEPTH_Q,
      F_STAIRS_HEIGHT,
      F_STAIRS_LENGTH,

      _PARAMS_COUNT
    };

    std::vector<real> wallM;
    std::vector<real> windowsM;
    std::vector<real> intM;
    std::vector<real> woodenM;
    std::vector<real> metalM;
    std::vector<real> roofM;
    std::vector<real> doorM;

    //functions to create sections. They are called with the same arguments and use too many of them,
    //so they use input generator parameters directly

    auto params_to_window_split_type = [&](real v) -> WindowSplit
    {
      if (v<0.5)
        return WindowSplit::WS_LEFT;
      if (v<1.5)
        return WindowSplit::WS_RIGHT;
      if (v<2.5)
        return WindowSplit::WS_NONE;
      else
        return WindowSplit::WS_BOTH;
    };

    auto section_code_to_types_vec = [](real code, real count) -> std::vector<real>
    {
      real count_t = 2;
      for (int i=0; i<count-1; i++)
        count_t *= 4;
      
      std::vector<real> tmp_res;
      while (code >= 2*count_t)
          code -= 2*count_t;
      while (count_t >= 1)
      {
        tmp_res.push_back(code >= count_t);
        if (code >= count_t)
          code -= count_t;
        count_t /= 2;
      }

      std::vector<real> res;
      for (int i=0;i<tmp_res.size();i+=2)
      {
        res.push_back(2*tmp_res[i] + tmp_res[i+1]);
      }

      return res;
    };

    auto calculate_building_width = [&](real section_code, real stripes_count, real wall_stripe_size) -> real
    {
      std::vector<real> stripe_types = section_code_to_types_vec(section_code, stripes_count);
      real sz = 0;
      for (int i=0; i<stripes_count; i++)
      {
        real stripe_type = stripe_types[i];
        if (stripe_type < 0.5)
          sz += wall_stripe_size;
        else if (stripe_type < 1.5)
          sz += params[F_BASE_WINDOW_STRIPE_SIZE];
        else if (stripe_type < 2.5)
          sz += params[F_BALCONY_WINDOW_STRIPE_SIZE];
      }
      return sz;
    };

    auto make_section = [&](bool left_end, bool right_end, vec3 start_point, real walls_dist, real height,
                            real section_code, real stripes_count, real wall_stripe_size,
                            bool is_side, vec3 vec_x, vec3 vec_z,
                            WindowQuality window_quality,
                            BalconyQuality balcony_quality) -> vec3
    {
      std::vector<real> stripe_types = section_code_to_types_vec(section_code, stripes_count);
      for (int i=0; i<stripes_count; i++)
      {
        real stripe_type = right_end ? stripe_types[stripe_types.size() - i - 1] : stripe_types[i];
        real stripe_size = 0;

        if (stripe_type < 0.5)
        {
          //wall stripe
          make_stripe_wall(wallM, intM, window_quality >= WindowQuality::WQ_HIGH,
                           params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                           Quad(start_point, wall_stripe_size*vec_x, vec3(0, height, 0), vec_z));
          make_stripe_wall(wallM, intM, window_quality >= WindowQuality::WQ_HIGH,
                           params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                           Quad(start_point - walls_dist*vec_z, wall_stripe_size*vec_x, vec3(0, height, 0), -vec_z));
          stripe_size = wall_stripe_size;
        }
        else if (stripe_type < 1.5)
        {
          //windows stripe
          make_stripe_windows(windowsM, wallM, woodenM, intM,
                  window_quality, params[I_BASE_WINDOW_HOR_SPLITS], params[I_BASE_WINDOW_VERT_SPLITS], 
                  params_to_window_split_type(params[I_BASE_WINDOW_SPLIT_TYPE]), params[F_BASE_WINDOW_SPLIT_H],
                  params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                  params[F_BASE_WINDOW_WIDTH_Q], params[F_BASE_WINDOW_HEIGHT_Q],
                  params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                  Quad(start_point, params[F_BASE_WINDOW_STRIPE_SIZE]*vec_x, vec3(0, height, 0), vec_z));
          make_stripe_windows(windowsM, wallM, woodenM, intM,
                  window_quality, params[I_BASE_WINDOW_HOR_SPLITS], params[I_BASE_WINDOW_VERT_SPLITS], 
                  params_to_window_split_type(params[I_BASE_WINDOW_SPLIT_TYPE]), params[F_BASE_WINDOW_SPLIT_H],
                  params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                  params[F_BASE_WINDOW_WIDTH_Q], params[F_BASE_WINDOW_HEIGHT_Q],
                  params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                  Quad(start_point - walls_dist*vec_z, params[F_BASE_WINDOW_STRIPE_SIZE]*vec_x, vec3(0, height, 0), -vec_z));
          stripe_size = params[F_BASE_WINDOW_STRIPE_SIZE];
        }
        else if (stripe_type < 2.5)
        {
          //balcony stripe
        make_stripe_balcony(windowsM, wallM, woodenM, intM, metalM,
                  window_quality, params[I_BALCONY_WINDOW_HOR_SPLITS], params[I_BALCONY_WINDOW_VERT_SPLITS], 
                  params_to_window_split_type(params[I_BALCONY_WINDOW_SPLIT_TYPE]), params[F_BALCONY_WINDOW_SPLIT_H],
                  params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                  params[F_BALCONY_WINDOW_WIDTH_Q], params[F_BALCONY_WINDOW_HEIGHT_Q],
                  params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q],
                  balcony_quality, 
                  params[I_FLOORS_COUNT], params[F_BALCONY_HEIGHT_Q], params[F_BALCONY_DEPTH_Q], 
                  params[F_BALCONY_BARS_THICKNESS], params[F_BALCONY_BARS_DISTANCE], params[I_BALCONY_START_FLOOR],
                  Quad(start_point, params[F_BALCONY_WINDOW_STRIPE_SIZE]*vec_x, vec3(0, height, 0), vec_z));
        make_stripe_balcony(windowsM, wallM, woodenM, intM, metalM,
                  window_quality, params[I_BALCONY_WINDOW_HOR_SPLITS], params[I_BALCONY_WINDOW_VERT_SPLITS], 
                  params_to_window_split_type(params[I_BALCONY_WINDOW_SPLIT_TYPE]), params[F_BALCONY_WINDOW_SPLIT_H],
                  params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                  params[F_BALCONY_WINDOW_WIDTH_Q], params[F_BALCONY_WINDOW_HEIGHT_Q],
                  params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], 
                  balcony_quality,
                  params[I_FLOORS_COUNT], params[F_BALCONY_HEIGHT_Q], params[F_BALCONY_DEPTH_Q], 
                  params[F_BALCONY_BARS_THICKNESS], params[F_BALCONY_BARS_DISTANCE], params[I_BALCONY_START_FLOOR],
                  Quad(start_point - walls_dist*vec_z, params[F_BALCONY_WINDOW_STRIPE_SIZE]*vec_x, vec3(0, height, 0), -vec_z));
        stripe_size = params[F_BALCONY_WINDOW_STRIPE_SIZE];
        }
        else
        {
          //we shouldn't be here
          assert(false);
        }

        if (!is_side)
        {
          if (window_quality >= WindowQuality::WQ_HIGH)
          {
            make_insides(intM, walls_dist, !(left_end && (i == 0)), params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                        Quad(start_point, vec3(stripe_size,0,0), vec3(0, height, 0), vec3(0, 0, 1)));
          }
          make_roof_segment(wallM, roofM, params[F_ROOF_BASE_H], PI/(3 + params[F_ROOF_SLOPE]), 
                            left_end && (i == 0), (right_end) && (i == stripe_types.size()-1), 
                            Quad(start_point + vec3(0, height, params[F_ROOF_OVERSIZE]), vec3(stripe_size,0,0), 
                                vec3(0,0,-walls_dist - 2*params[F_ROOF_OVERSIZE]), vec3(0,1,0)));
        }
        start_point += vec_x*stripe_size;
      }
      return start_point;
    };

    auto make_entrance_section = [&](vec3 start_point, real walls_dist, real height,
                                     WindowQuality window_quality) -> vec3
    {
      make_stripe_entrance(windowsM, wallM, woodenM, intM, metalM, 
                          window_quality, params[I_ENTRANCE_WINDOW_HOR_SPLITS], params[I_ENTRANCE_WINDOW_VERT_SPLITS], 
                          params_to_window_split_type(params[I_ENTRANCE_WINDOW_SPLIT_TYPE]), params[F_ENTRANCE_WINDOW_SPLIT_H],
                          params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                          params[F_ENTRANCE_WINDOW_WIDTH_Q], params[F_ENTRANCE_WINDOW_HEIGHT_Q],
                          params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                          params[F_DOOR_OFFSET_Q], params[F_DOOR_WIDTH_Q], params[F_DOOR_HEIGHT_Q],
                          params[F_DOOR_DEPTH_Q], params[F_STAIRS_HEIGHT], params[F_STAIRS_LENGTH],
                          Quad(start_point, vec3(params[F_ENTRANCE_WINDOW_STRIPE_SIZE],0,0), vec3(0, height, 0), vec3(0, 0, 1)));
      make_stripe_windows(windowsM, wallM, woodenM, intM,
                          window_quality, params[I_ENTRANCE_WINDOW_HOR_SPLITS], params[I_ENTRANCE_WINDOW_VERT_SPLITS], 
                          params_to_window_split_type(params[I_ENTRANCE_WINDOW_SPLIT_TYPE]), params[F_ENTRANCE_WINDOW_SPLIT_H],
                          params[F_WINDOW_OUTER_FRAME_WIDTH], params[F_WINDOW_INNER_FRAME_WIDTH], params[F_WINDOW_GLASS_DEPTH],
                          params[F_ENTRANCE_WINDOW_WIDTH_Q], params[F_ENTRANCE_WINDOW_HEIGHT_Q],
                          params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                          Quad(start_point - vec3(0,0,walls_dist), vec3(params[F_ENTRANCE_WINDOW_STRIPE_SIZE],0,0), vec3(0, height, 0), vec3(0, 0, -1)));
      real stripe_size = params[F_ENTRANCE_WINDOW_STRIPE_SIZE];

      if (window_quality >= WindowQuality::WQ_HIGH)
      {
        make_insides(intM, walls_dist, true, params[F_WALL_THICKNESS], params[F_BOTTOM_OFFSET_Q], params[F_TOP_OFFSET_Q], params[I_FLOORS_COUNT],
                    Quad(start_point, vec3(stripe_size,0,0), vec3(0, height, 0), vec3(0, 0, 1))); 
      }
      make_roof_segment(wallM, roofM, params[F_ROOF_BASE_H], PI/(3 + params[F_ROOF_SLOPE]), 
                        false, false, 
                        Quad(start_point + vec3(0, height, params[F_ROOF_OVERSIZE]), vec3(stripe_size,0,0), 
                             vec3(0,0,-walls_dist - 2*params[F_ROOF_OVERSIZE]), vec3(0,1,0)));

      start_point.x += stripe_size;
      return start_point;
    };

    for (int i=0;i<params.size();i++)
    {
      std::cerr<<params[i]<<", ";
    }
    std::cerr<<"\n";

    WindowQuality wq = quality.quality_level >= ModelQuality::HIGH ? WQ_HIGH : WQ_LOW;
    BalconyQuality bq = quality.quality_level >= ModelQuality::HIGH ? BQ_HIGH : (quality.quality_level >= ModelQuality::MEDIUM ? BQ_LOW : BQ_NONE);

    real size_w = calculate_building_width(params[I_SS_CODE], params[I_SS_STRIPES], params[F_SS_WALL_STRIPE_SIZE]);
    real height = (params[I_FLOORS_COUNT] + params[F_BOTTOM_OFFSET_Q] + params[F_TOP_OFFSET_Q]);
    real start_x = 0;
    vec3 start_point = vec3(start_x,0,0);
    start_point = make_section(true, false, start_point, size_w, height, params[I_ES_CODE], params[I_ES_STRIPES], params[F_ES_WALL_STRIPE_SIZE],
                               false, vec3(1,0,0), vec3(0,0,1), wq, bq);
    for (int i=0;i<params[I_ENTRANCES_COUNT] - 1;i++)
    {
      start_point = make_entrance_section(start_point, size_w, height, wq);
      start_point = make_section(false, false, start_point, size_w, height, params[I_MS_CODE], params[I_MS_STRIPES], params[F_MS_WALL_STRIPE_SIZE],
                                 false, vec3(1,0,0), vec3(0,0,1), wq, bq);
    }
    start_point = make_entrance_section(start_point, size_w, height, wq);
    start_point = make_section(false, true, start_point, size_w, height, params[I_ES_CODE], params[I_ES_STRIPES], params[F_ES_WALL_STRIPE_SIZE],
                               false, vec3(1,0,0), vec3(0,0,1), wq, bq);

    real size_l = start_point.x - start_x;
    make_section(false, false, vec3(start_x,0,-size_w), size_l, height, params[I_SS_CODE], params[I_SS_STRIPES], params[F_SS_WALL_STRIPE_SIZE],
                               true, vec3(0,0,1), vec3(-1,0,0), wq, bq);

    std::vector<std::vector<real> *> models = {&wallM, &windowsM, &intM, &woodenM, &metalM, &roofM, &doorM};
    std::vector<std::string> names = {"main_part", "windows", "interior", "wooden_parts", "metal_parts", "roof", "door"};
    int sz = 0;
    for (int i=0;i<models.size(); i++)
      sz += models[i]->size();
    model.resize(sz);

    int k = 0;
    PartOffsets po;
    for (int i=0;i<models.size(); i++)
    {
      debug("%s: %d vertices\n", names[i].c_str(), models[i]->size());
      if (!(models[i]->empty()))
        po.push_back({names[i], k});
      for (real &v : *(models[i]))
      {
        model[k] = v;
        k++;
      }
    }

    real t = height > size_l ? 1/height : 1/size_l;
    dmat43 sc2 = scale(ident<dfloat>(), dvec3{t, t, t});
    transform(model, sc2);

    return po;
  }
    
  PartOffsets create_building_2(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality)
  {
    return create_building_2t<dfloat>(params, model, quality);
  }
}
