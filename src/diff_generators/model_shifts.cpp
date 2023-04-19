#pragma once
#include <vector>
#include "vectors.h"
#include <string>
#include "diff_geometry_generation.h"

namespace dgen
{
  std::vector<float> shift_model(std::vector<float> model, float x, float y, float z)
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      model[i * FLOAT_PER_VERTEX] += x;
      model[i * FLOAT_PER_VERTEX + 1] += y;
      model[i * FLOAT_PER_VERTEX + 2] += z;
    }
    return model;
  }

  std::vector<float> matrix_mul_model(std::vector<float> model, float A[3][3])
  {
    for (int i = 0; i < model.size() / FLOAT_PER_VERTEX; ++i)
    {
      float x = 0, y = 0, z = 0, a = 0, b = 0, c = 0;
      for (int j = 0; j < 3; ++j)
      {
        x += A[0][j] * model[i * FLOAT_PER_VERTEX + j];
        y += A[1][j] * model[i * FLOAT_PER_VERTEX + j];
        z += A[2][j] * model[i * FLOAT_PER_VERTEX + j];
        a += A[0][j] * model[i * FLOAT_PER_VERTEX + j + 3];
        b += A[1][j] * model[i * FLOAT_PER_VERTEX + j + 3];
        c += A[2][j] * model[i * FLOAT_PER_VERTEX + j + 3];
      }
      model[i * FLOAT_PER_VERTEX] = x;
      model[i * FLOAT_PER_VERTEX + 1] = y;
      model[i * FLOAT_PER_VERTEX + 2] = z;
      model[i * FLOAT_PER_VERTEX + 3] = a;
      model[i * FLOAT_PER_VERTEX + 4] = b;
      model[i * FLOAT_PER_VERTEX + 5] = c;
    }
    return model;
  }
}