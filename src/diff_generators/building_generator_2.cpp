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
                                 WindowSplit ws, real split_h_q, real thick,
                                 real outer_frame_width, real inner_frame_width,
                                 real glass_deep_q,
                                 real window_width_q, real window_height_q,
                                 real bottom_offset_q, real top_offset_q, real floors_count,
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
                                    WindowSplit ws, real split_h_q, real thick,
                                    real outer_frame_width, real inner_frame_width,
                                    real glass_deep_q,
                                    real window_width_q, real window_height_q,
                                    real bottom_offset_q, real top_offset_q, real floors_count,
                                    real door_offset_q, real door_width_q, real door_height_q, real door_deep_q,
                                    real stairs_h_q, real stairs_w_q,
                                    const Quad &q)
    {
      bool has_interior = (wq == WindowQuality::WQ_HIGH);
      real y_mul = 1/(floors_count + bottom_offset_q + top_offset_q);
      Quad q_part = q;

      //door
      
      real w_off = 0.5f * (1 - door_width_q) * length(y_mul*q.v2) / length(q.v1);
      real door_wq = door_width_q * length(y_mul*q.v2) / length(q.v1);
      Quad q_door(q.p1 + door_offset_q*y_mul*q.v2 + w_off*q.v1 - door_deep_q*thick*q.n, vec2(0,0), door_wq * q.v1, vec2(1,0), door_height_q*y_mul*q.v2, vec2(0,1),q.n);
      make_door(M_door, 0.2*thick, q_door);
      
      //doorway
      {
        Quad q_l(q.p1 + w_off*q.v1, vec2(0,0), -thick*q.n, vec2(1,0), (bottom_offset_q+1) * y_mul * q.v2, vec2(0, 1), q.v1);
        Quad q_r(q.p1 + (w_off + door_wq)*q.v1, vec2(0,0), (bottom_offset_q+1) * y_mul * q.v2, vec2(0, 1), -thick*q.n, vec2(1,0), -q.v1);
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
      
        Quad q_doorway_above(q.p1 + (door_offset_q+door_height_q)*y_mul*q.v2 + w_off*q.v1, vec2(0,0), -thick*q.n, vec2(1,0), door_wq * q.v1, vec2(0, 1), -q.v2);
        make_quad(M_int, q_doorway_above);
      }

      //box under the door and stairs
      {
        Quad q_under(q.p1 + w_off*q.v1, vec2(0,0), door_wq * q.v1, vec2(1,0), door_offset_q * y_mul * q.v2, vec2(0, 1), q.n);
        make_box(M_int, thick, q_under);

        real sw = stairs_w_q * length(y_mul*q.v2);
        Quad s_quad(q.p1 + w_off*q.v1, vec2(0,0), door_wq * q.v1, vec2(1,0), q.v2, vec2(0, 1), q.n);
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

    auto make_balcony = [&](std::vector<real> &M_plate, std::vector<real> &M_bars, 
                            real plate_thick, real depth, real bars_thick, real bars_distance, real bars_height,
                            const Quad &q)
    {
      Quad q_plate(q.p1 + depth*q.n, vec2(0,0), q.v1, vec2(1,0), -depth*q.n, vec2(0,1), q.v2);
      make_box(M_plate, plate_thick, q_plate);

      real t = 0;
      real b_th = bars_thick*length(q.v1);
      real b_th_y = bars_thick*length(q.v1) / length(q.v2);
      //bars in front
      while (t < 1)
      {
        make_box(M_bars, b_th, Quad(q.p1 + depth*q.n + t*q.v1, vec2(0,0), bars_thick*q.v1, vec2(1,0), bars_height*q.v2, vec2(0,1), q.n));
        t += bars_distance;
      }
      make_box(M_bars, b_th, Quad(q.p1 + depth*q.n + bars_height*q.v2, vec2(0,0), q.v1, vec2(1,0), b_th_y*q.v2, vec2(0,1), q.n));

      t = 0;
      real b_dist = bars_distance*length(q.v1);
      while (t < depth)
      {
        make_box(M_bars, b_th, Quad(q.p1 + (depth-t)*q.n, vec2(0,0), bars_thick*q.v1, vec2(1,0), bars_height*q.v2, vec2(0,1), q.n));
        make_box(M_bars, b_th, Quad(q.p1 + (depth-t)*q.n + (1-bars_thick)*q.v1, vec2(0,0), bars_thick*q.v1, vec2(1,0), bars_height*q.v2, vec2(0,1), q.n));
        t+= b_dist;
      }
      make_box(M_bars, depth - b_th, Quad(q.p1 + (depth-b_th)*q.n + bars_height*q.v2, vec2(0,0), bars_thick*q.v1, vec2(1,0), b_th_y*q.v2, vec2(0,1), q.n));
      make_box(M_bars, depth - b_th, Quad(q.p1 + (depth-b_th)*q.n + (1-bars_thick)*q.v1 + bars_height*q.v2, vec2(0,0), bars_thick*q.v1, vec2(1,0), b_th_y*q.v2, vec2(0,1), q.n));
    };

    auto make_stripe_balcony =  [&](std::vector<real> &M_glass, std::vector<real> &M_wall,
                                    std::vector<real> &M_frame, std::vector<real> &M_int,
                                    std::vector<real> &M_bars,
                                    WindowQuality wq, real hor_parts, real vert_parts,
                                    WindowSplit ws, real split_h_q, real thick,
                                    real outer_frame_width, real inner_frame_width,
                                    real glass_deep_q,
                                    real window_width_q, real window_height_q,
                                    real bottom_offset_q, real top_offset_q, real floors_count,
                                    real balcony_height_q, real balcony_depth_q,
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
          make_balcony(M_int, M_bars, 0.04 * length(y_mul*q.v2), balcony_depth_q * length(y_mul*q.v2), bars_thick, bars_distance, balcony_height_q, q_part);

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

    std::vector<real> wallM;
    std::vector<real> windowsM;
    std::vector<real> intM;
    std::vector<real> woodenM;
    std::vector<real> metalM;

    make_stripe_entrance(windowsM, wallM, woodenM, intM, metalM, 
                        WindowQuality::WQ_HIGH, 1, 2, 
                        WindowSplit::WS_BOTH, 0.6, 0.1,
                        0.1, 0.07, 0.1,
                        0.7, 0.7,
                        0.8, 0.2, 4,
                        0.6, 0.6, 1, 0.15,
                        0.1, 0.13,
                        Quad(vec3(0.0, 0.0, 0), vec2(0,0), vec3(0.2,0,0), vec2(1,0), vec3(0, 1, 0), vec2(0, 1), vec3(0, 0, 1)));
    
    make_stripe_balcony(windowsM, wallM, woodenM, intM, metalM,
                  WindowQuality::WQ_HIGH, 1, 3, 
                  WindowSplit::WS_LEFT, 0.6, 0.1,
                  0.1, 0.07, 0.1,
                  0.7, 0.7,
                  0.8, 0.2, 4, 
                  0.4, 0.5, 0.015, 0.05, 1,
                  Quad(vec3(0.2, 0.0, 0), vec2(0,0), vec3(0.3,0,0), vec2(1,0), vec3(0, 1, 0), vec2(0, 1), vec3(0, 0, 1)));
    for (int i=0;i<2;i++)
    {
      make_stripe_windows(windowsM, wallM, woodenM, intM,
                  WindowQuality::WQ_HIGH, 1, 2, 
                  WindowSplit::WS_BOTH, 0.6, 0.1,
                  0.1, 0.07, 0.1,
                  0.7, 0.7,
                  0.8, 0.2, 4,
                  Quad(vec3(0.2*i + 0.5, 0.0, 0), vec2(0,0), vec3(0.2,0,0), vec2(1,0), vec3(0, 1, 0), vec2(0, 1), vec3(0, 0, 1)));
    }
    make_stripe_wall(wallM, intM, true, 0.1, 0.8, 0.2, 4,
                     Quad(vec3(0.9, 0.0, 0), vec2(0,0), vec3(0.2,0,0), vec2(1,0), vec3(0, 1, 0), vec2(0, 1), vec3(0, 0, 1)));


    std::vector<std::vector<real> *> models = {&wallM, &windowsM, &intM, &woodenM, &metalM};
    
    int sz = 0;
    for (int i=0;i<models.size(); i++)
      sz += models[i]->size();
    model.resize(sz);

    int k = 0;
    for (int i=0;i<models.size(); i++)
    {
      for (real &v : *(models[i]))
      {
        model[k] = v;
        k++;
      }
    }

    return {
        {"main_part", 0},
        {"windows", wallM.size()},
        {"interior", windowsM.size() + wallM.size()},
        {"wooden_parts", windowsM.size() + wallM.size() + intM.size()},
        {"metal_parts", windowsM.size() + wallM.size() + intM.size() + woodenM.size()}};
  }
    
  PartOffsets create_building_2(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality)
  {
    return create_building_2t<dfloat>(params, model, quality);
  }
}
