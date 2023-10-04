#include "render_target.h"
#include "engine.h"

bool RenderTarget::create(int w, int h)
{
    width = w;
    height = h;
    float borderColor[] = {0.0f, 0.0f, 0.0f, 0.0f};
    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    frBuffer = create_framebuffer();

    tex = engine::textureManager->create_texture(width, height, texFmt, 1);
    glBindTexture(GL_TEXTURE_2D, tex.texture);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex.texture, 0);

    unsigned int attachments[1] = { GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
        return false;
    }
    else
    {
        debugl(10, "Render target created %d %d", width, height);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
    return true;
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
  engine::textureManager->delete_tex(tex);
  delete_framebuffer(frBuffer);
}