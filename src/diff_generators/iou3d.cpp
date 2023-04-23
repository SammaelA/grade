#include "iou3d.h"
#include <cstdlib>
#include <ctime>

bool check_tr_line(double x1, double y1, double z1,
  double x2, double y2, double z2,
  double x3, double y3, double z3,
  double x, double y, double z,
  double x_v, double y_v, double z_v)
{
  double e1[3] = {x2 - x1, y2 - y1, z2 - z1};
  double e2[3] = {x3 - x1, y3 - y1, z3 - z1};
  double T[3] = {x - x1, y - y1, z - z1};
  double P[3] = {y_v * e2[2] - z_v * e2[1], z_v * e2[0] - x_v * e2[2], x_v * e2[1] - y_v * e2[0]};
  double Q[3] = {T[1] * e1[2] - T[2] * e1[1], T[2] * e1[0] - T[0] * e1[2], T[0] * e1[1] - T[1] * e1[0]};
  double dots[4] = {Q[0] * e2[0] + Q[1] * e2[1] + Q[2] * e2[2], P[0] * T[0] + P[1] * T[1] + P[2] * T[2], Q[0] * x_v + Q[1] * y_v + Q[2] * z_v, P[0] * e1[0] + P[1] * e1[1] + P[2] * e1[2]};
  if (dots[3] == 0)
  {
    return false;
  }
  return dots[0] / dots[3] >= 0 && dots[1] / dots[3] >= 0 && dots[2] / dots[3] >= 0 && dots[1] / dots[3] + dots[2] / dots[3] <= 1;
}

/*bool inter_tr_line(double x1, double y1, double z1,
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
}*/

double iou3d(const std::vector<float> &model1, const std::vector<float> &model2, 
            double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, double step)
{
  double intersect = 0, interconnect = 0;
  double result = 0;
  /*int X = (x_e - x_s) / step + 1;
  int Y = (y_e - y_s) / step + 1;
  int Z = (z_e - z_s) / step + 1;
  glm::vec4 *arr = new glm::vec4[X * Y * Z];
  int a = 0;*/
  for (double x = x_s; x <= x_e; x += step)
  {
    //int b = 0;
    for (double y = y_s; y <= y_e; y += step)
    {
      //int c = 0;
      for (double z = z_s; z <= z_e; z += step)
      {
        int num1 = 0, num2 = 0;
        double x_vec = rand() % 2 + 0.5;
        double y_vec = rand() % 2 + 0.5;
        double z_vec = rand() % 2 + 0.5;
        for (int i = 0; i < model1.size(); i += 3 * FLOAT_PER_VERTEX)
        {
          if (check_tr_line(model1[i], model1[i + 1], model1[i + 2], 
              model1[i + FLOAT_PER_VERTEX], model1[i + 1 + FLOAT_PER_VERTEX], model1[i + 2 + FLOAT_PER_VERTEX], 
              model1[i + FLOAT_PER_VERTEX * 2], model1[i + 1 + FLOAT_PER_VERTEX * 2], model1[i + 2 + FLOAT_PER_VERTEX * 2], 
              x, y, z, x_vec, y_vec, z_vec))
          {
            ++num1;
          }
        }
        for (int i = 0; i < model2.size(); i += 3 * FLOAT_PER_VERTEX)
        {
          if (check_tr_line(model2[i], model2[i + 1], model2[i + 2], 
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
          //arr[X * Y * c + X * b + a] = glm::vec4{1, 0, 0, 0.5};
        }/* else if (prod1 || prod2)
        {
          arr[X * Y * c + X * b + a] = glm::vec4{0, 1, 1, 0.5};
        }
        else
        {
          arr[X * Y * c + X * b + a] = glm::vec4{0, 0, 0, 0};
        }*/
        if (prod1 || prod2)
        {
          interconnect += 1;
        }
        //++c;
      }
      //++b;
    }
    //++a;
  }
  /*glm::vec3 p0 = {x_s, y_s, z_s};
  glm::vec3 p1 = {x_e, y_e, z_e};
  glm::ivec3 vox_count = {X, Y, Z};
  VoxelArray <glm::vec4> vox(p0, p1, vox_count, {0, 0, 0, 0}, arr);
  voxelization::render_test_3d(vox);*/
  if (interconnect == 0)
  {
    return 0;
  }
  return intersect / interconnect;
}
