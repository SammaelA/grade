#include "graphics_utils.h"
#include "common_utils/utility.h"

void print_FB_status(GLuint status)
{
    if (status == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT\n");
    else if (status == GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS\n");
    else if (status == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT\n");
    else if (status == GL_FRAMEBUFFER_UNSUPPORTED)
        debugl(9,"GL_FRAMEBUFFER_UNSUPPORTED\n");
    else  debugl(9,"GL_FRAMEBUFFER_INCOMPLETE %#010x\n",status);
}