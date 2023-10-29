#include "custom_diff_render/custom_diff_render.h"
#include "common_utils/utility.h"

IDiffRender *get_halfgpu_custom_diff_render()
{
  logerr("custom_diff_render module is not included. Change settings in CMakeLists.txt to use it.");
  return nullptr;
}