#include "tables_generator.h"
#include "diff_geometry_generation.h"
#include <cppad/cppad.hpp>

namespace dgen
{
  inline std::vector<float> get_cube()
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
      1.00,0.00,0.00, 0.00,-1.00,0.00, 1.00,0.00,
      0.00,0.00,1.00, 0.00,-1.00,0.00, 0.00,0.00,
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
  void create_simple_table(const std::vector<dfloat> &params, std::vector<dfloat> &model, ModelQuality quality)
  {
    dfloat countertop_w = params[0];
    dfloat countertop_h = params[1];
    dfloat countertop_thick = params[2];
    dfloat leg_h = params[3];
    dfloat leg_thick = params[4];
    dfloat leg_offset = params[5];

    {
      std::vector<float> m1 = get_cube();
      std::vector<dfloat> m;
      for (auto &v : m1)
        m.push_back(v);
      dmat43 mat = mul(translate(ident(), dvec3{-0.5*countertop_w, leg_h, -0.5*countertop_h}),
                       scale(ident(), dvec3{countertop_w, countertop_thick, countertop_h}));
      transform(m, mat);
      for (auto &v : m)
        model.push_back(v);
    }

    {
      std::vector<float> m1 = get_cube();
      std::vector<dfloat> m;
      for (auto &v : m1)
        m.push_back(v);
      dmat43 mat = mul(translate(ident(), dvec3{-(0.5*countertop_w - leg_offset), 0, -(0.5*countertop_h - leg_offset)}),
                       scale(ident(), dvec3{leg_thick, leg_h, leg_thick}));
      transform(m, mat);
      for (auto &v : m)
        model.push_back(v);
    }

    {
      std::vector<float> m1 = get_cube();
      std::vector<dfloat> m;
      for (auto &v : m1)
        m.push_back(v);
      dmat43 mat = mul(translate(ident(), dvec3{(0.5*countertop_w - leg_offset - leg_thick), 0, -(0.5*countertop_h - leg_offset)}),
                       scale(ident(), dvec3{leg_thick, leg_h, leg_thick}));
      transform(m, mat);
      for (auto &v : m)
        model.push_back(v);
    }

    {
      std::vector<float> m1 = get_cube();
      std::vector<dfloat> m;
      for (auto &v : m1)
        m.push_back(v);
      dmat43 mat = mul(translate(ident(), dvec3{-(0.5*countertop_w - leg_offset), 0, (0.5*countertop_h - leg_offset - leg_thick)}),
                       scale(ident(), dvec3{leg_thick, leg_h, leg_thick}));
      transform(m, mat);
      for (auto &v : m)
        model.push_back(v);
    }

    {
      std::vector<float> m1 = get_cube();
      std::vector<dfloat> m;
      for (auto &v : m1)
        m.push_back(v);
      dmat43 mat = mul(translate(ident(), dvec3{(0.5*countertop_w - leg_offset - leg_thick), 0, (0.5*countertop_h - leg_offset - leg_thick)}),
                       scale(ident(), dvec3{leg_thick, leg_h, leg_thick}));
      transform(m, mat);
      for (auto &v : m)
        model.push_back(v);
    }
  }
}