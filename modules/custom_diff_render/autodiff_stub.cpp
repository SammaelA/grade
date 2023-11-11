#include "autodiff.h"
#include "common_utils/utility.h"

void __enzyme_fwddiff(void*, ...)
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
}
int enzyme_dupnoneed = 0;
int enzyme_dup = 0;
float __enzyme_v[1024] = {0};