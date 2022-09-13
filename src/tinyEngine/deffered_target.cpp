#include "deffered_target.h"
#include "tinyEngine/engine.h"

bool DefferedTarget::create(int w, int h)
{
    width = w;
    height = h;
    float borderColor[] = {0.0f, 0.0f, 0.0f, 0.0f};
    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    frBuffer = create_framebuffer();

    depthTex = engine::textureManager->create_texture(width, height, GL_DEPTH_COMPONENT16, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);
    glBindTexture(GL_TEXTURE_2D, depthTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);

    colorTex = engine::textureManager->create_texture(width, height, colorTexFmt, 1, NULL, GL_RGB, GL_UNSIGNED_BYTE);
    glBindTexture(GL_TEXTURE_2D, colorTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    normalsTex = engine::textureManager->create_texture(width, height, normalsTexFmt, 1, NULL, GL_RGB, GL_UNSIGNED_BYTE);
    glBindTexture(GL_TEXTURE_2D, normalsTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    viewPosTex = engine::textureManager->create_texture(width, height, viewPosTexFmt, 1, NULL, GL_RGB, GL_UNSIGNED_BYTE);
    glBindTexture(GL_TEXTURE_2D, viewPosTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    worldPosTex = engine::textureManager->create_texture(width, height, worldPosTexFmt, 1, NULL, GL_RGB, GL_UNSIGNED_BYTE);
    glBindTexture(GL_TEXTURE_2D, worldPosTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex.texture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorTex.texture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, normalsTex.texture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, viewPosTex.texture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, worldPosTex.texture, 0);

    unsigned int attachments[4] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3};
    glDrawBuffers(4, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    else
    {
        debugl(10, "Deferred target created %d %d", width, height);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
}
void DefferedTarget::target()
{
  glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
  glViewport(0, 0, width, height);
  glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
DefferedTarget::~DefferedTarget()
{
  engine::textureManager->delete_tex(colorTex);
  engine::textureManager->delete_tex(depthTex);
  engine::textureManager->delete_tex(normalsTex);
  engine::textureManager->delete_tex(viewPosTex);
  engine::textureManager->delete_tex(worldPosTex);
  delete_framebuffer(frBuffer);
}