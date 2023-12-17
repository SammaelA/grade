#include "autodiff.h"
#include "common_utils/utility.h"

void __enzyme_fwddiff(void*, ...)
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
}
float __enzyme_autodiffFloat(float (*)(float), float)
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
}

double __enzyme_autodiffDouble(double (*)(double), double)
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
}

float __enzyme_autodiff(...)
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
}

int enzyme_dupnoneed = 0;
int enzyme_dup = 0;
float __enzyme_v[1024] = {0};