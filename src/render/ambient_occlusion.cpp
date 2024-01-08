#include "ambient_occlusion.h"
static glm::vec2 calcFocalLen(float fovRad, float width, float height)
{
    glm::vec2 res;
    res.x = 1.f/tanf(0.5 * fovRad)*(height/width);
    res.y = 1.f/tanf(0.5 * fovRad);
    return res;
}
HBAORenderer::HBAORenderer():
              noise(engine::textureManager->get("colored_noise")),
              shader("hbao.fs")
{

}
HBAORenderer::~HBAORenderer()
{
  engine::textureManager->delete_tex(aoTex);
  delete_framebuffer(frBuffer);
}
void HBAORenderer::render(float fovRad, float width, float height, GLuint viewPosTex)
{
    glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
    glViewport(0, 0, width, height);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    shader.use();
    shader.get_shader().texture("noiseTex",noise.texture);
    shader.get_shader().texture("viewPosTex",viewPosTex);
    shader.get_shader().uniform("FocalLen",calcFocalLen(fovRad,width,height));
    shader.render();
    
    glEnable(GL_DEPTH_TEST);
}
void HBAORenderer::create(int w, int h)
{
    width = w;
    height = h;
    float borderColor[] = {0.0f, 0.0f, 0.0f, 0.0f};
    frBuffer = create_framebuffer();

    aoTex = engine::textureManager->create_texture(width, height, aoTexFmt, 1, nullptr, GL_RGB, GL_UNSIGNED_BYTE);
    glBindTexture(GL_TEXTURE_2D, aoTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, frBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, aoTex.texture, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}
