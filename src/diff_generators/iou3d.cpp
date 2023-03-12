#include "iou3d.h"
#include <cstdlib>
#include <ctime>

bool inter_tr_line(double x1, double y1, double z1,
  double x2, double y2, double z2,
  double x3, double y3, double z3,
  double x, double y, double z,
  double x_v, double y_v, double z_v)
{
  double matr[4][3];
  matr[0][0] = x1 - x2;
  matr[0][1] = y1 - y2;
  matr[0][2] = z1 - z2;
  matr[1][0] = x1 - x3;
  matr[1][1] = y1 - y3;
  matr[1][2] = z1 - z3;
  matr[2][0] = x_v;
  matr[2][1] = y_v;
  matr[2][2] = z_v;
  matr[3][0] = x1 - x;
  matr[3][1] = y1 - y;
  matr[3][2] = z1 - z;
  double a, b, c;
  for (int k = 0; k < 3; ++k)
  {
    for (int i = k; i < 3; ++i)
    {
      if (matr[k][i] != 0)
      {
        double div = matr[k][i];
        for (int j = 0; j < 4; ++j)
        {
          double tmp = matr[j][i] / div;
          matr[j][i] = matr[j][k];
          matr[j][k] = tmp;
        }
        for (int j = 0; j < 3; ++j)
        {
          double div = matr[k][j];
          for (int u = 0; u < 4; ++u)
          {
            if (j != k && div != 0)
            {
              matr[u][j] = matr[u][j] / div - matr[u][k];
            }
          }
        }
        break;
      }
      if (i == 2)
      {
        return false;
      }
    }
  }
  return matr[3][0] >= 0 && matr[3][0] <= 1 && matr[3][1] >= 0 && matr[3][1] <= 1 && matr[3][1] + matr[3][0] <= 1 && matr[3][2] >= 0;
}

double iou3d(const std::vector<float> &model1, const std::vector<float> &model2, 
            double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, double step)
{
  double intersect = 0, interconnect = 0;
  double result = 0;
  for (double x = x_s; x <= x_e; x += step)
  {
    for (double y = y_s; y <= y_e; y += step)
    {
      for (double z = z_s; z <= z_e; z += step)
      {
        int num1 = 0, num2 = 0;
        double x_vec = rand() % 2 + 0.5;
        double y_vec = rand() % 2 + 0.5;
        double z_vec = rand() % 2 + 0.5;
        for (int i = 0; i < model1.size(); i += 3 * FLOAT_PER_VERTEX)
        {
          if (inter_tr_line(model1[i], model1[i + 1], model1[i + 2], 
              model1[i + FLOAT_PER_VERTEX], model1[i + 1 + FLOAT_PER_VERTEX], model1[i + 2 + FLOAT_PER_VERTEX], 
              model1[i + FLOAT_PER_VERTEX * 2], model1[i + 1 + FLOAT_PER_VERTEX * 2], model1[i + 2 + FLOAT_PER_VERTEX * 2], 
              x, y, z, x_vec, y_vec, z_vec))
          {
            ++num1;
          }
        }
        for (int i = 0; i < model2.size(); i += 3 * FLOAT_PER_VERTEX)
        {
          if (inter_tr_line(model2[i], model2[i + 1], model2[i + 2], 
              model2[i + FLOAT_PER_VERTEX], model2[i + 1 + FLOAT_PER_VERTEX], model2[i + 2 + FLOAT_PER_VERTEX], 
              model2[i + FLOAT_PER_VERTEX * 2], model2[i + 1 + FLOAT_PER_VERTEX * 2], model2[i + 2 + FLOAT_PER_VERTEX * 2], 
              x, y, z, x_vec, y_vec, z_vec))
          {
            ++num2;
          }
        }
        bool prod1 = (num1 % 2 == 1);
        bool prod2 = (num2 % 2 == 1);
        if (prod1 && prod2)
        {
          intersect += 1;
        }
        if (prod1 || prod2)
        {
          interconnect += 1;
        }
      }
    }
  }
  if (interconnect == 0)
  {
    return 0;
  }
  return intersect / interconnect;
}
