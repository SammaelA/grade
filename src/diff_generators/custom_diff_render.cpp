#include "custom_diff_render.h"
#include "common_utils/utility.h"

#ifdef USE_CUSTOM_DIFF_RENDER
#else
int custom_diff_render_main(int argc, char *argv[])
{
  logerr("custom differentiable renderer is not implemented!");
  return 0;
}
#endif