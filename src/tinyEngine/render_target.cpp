#include "render_target.h"

bool RenderTarget::create(int w, int h)
{
    width = w;
    height = h;
    float borderColor[] = {0.0f, 0.0f, 0.0f, 0.0f};
    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glGenFramebuffers(1, &frBuffer);

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, texFmt, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);

    unsigned int attachments[1] = { GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    else
    {
        debugl(10, "Render target created %d %d", width, height);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void RenderTarget::target()
{
  glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
  glViewport(0, 0, width, height);
  glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);
}
RenderTarget::~RenderTarget()
{
    glDeleteTextures(1, &tex);
    glDeleteFramebuffers(1, &frBuffer);
}