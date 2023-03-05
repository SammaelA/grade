#include "iou3d.h"
#include "diff_geometry_generation.h"
#include "common_utils/interpolation.h"
#include <cppad/cppad.hpp>

double iou3d(const std::vector<dfloat> &model1, const std::vector<dfloat> &model2, 
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
        double dist = -1;
        int num1 = 0, num2 = 0;
        double tmp = 0;
        for (int i = 0; i < model1.size(); i += FLOAT_PER_VERTEX)
        {
          tmp = pow(model1[i + 0] - x, 2) + pow(model1[i + 1] - y, 2) + pow(model1[i + 2] - z, 2);
          if (dist == -1 || tmp < dist || 
              (tmp == dist && (x - model1[i + 0]) * model1[i + 3] + (y - model1[i + 1]) * model1[i + 4] + (x - model1[i + 2]) * model1[i + 5] != 0))
          {
            dist = tmp;
            num1 = i;
          }
        }
        dist = -1;
        for (int i = 0; i < model2.size(); i += FLOAT_PER_VERTEX)
        {
          tmp = pow(model2[i + 0] - x, 2) + pow(model2[i + 1] - y, 2) + pow(model2[i + 2] - z, 2);
          if (dist == -1 || tmp < dist || 
              (tmp == dist && (x - model2[i + 0]) * model2[i + 3] + (y - model2[i + 1]) * model2[i + 4] + (x - model2[i + 2]) * model2[i + 5] != 0))
          {
            dist = tmp;
            num2 = i;
          }
        }
        bool prod1 = (x - model1[num1 + 0]) * model1[num1 + 3] + (y - model1[num1 + 1]) * model1[num1 + 4] + (x - model1[num1 + 2]) * model1[num1 + 5] <= 0;
        bool prod2 = (x - model2[num2 + 0]) * model2[num2 + 3] + (y - model2[num2 + 1]) * model2[num2 + 4] + (x - model2[num2 + 2]) * model2[num2 + 5] <= 0;
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