#pragma once
#include "differentiable_generators.h"


double iou3d(const std::vector<float> &model1, const std::vector<float> &model2, 
            double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, double step);
