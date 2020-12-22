#include "texture_atlas.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "tinyEngine/utility/model.h"
TextureAtlas::TextureAtlas(int w, int h):
colorTex(w,h),
mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{
    width = w;
    height = h;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind();   
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout<<"Framebuffer Incomplete"<<std::endl;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void TextureAtlas::set_grid(int w, int h)
{
    curNum = 0;
    isGrid = true;
    gridWN = width / w;
    gridHN = height / h;
}
int TextureAtlas::add_tex()
{
    curNum++;
    if (curNum>gridWN*gridHN)
        return -1;
    else
        return curNum - 1;
}
void TextureAtlas::process_tc_array(int num, std::vector<float> &tc)
{
    glm::vec2 tsc(1.0/gridWN,1.0/gridHN);
    glm::vec2 tsh(num % gridWN, num / gridWN);
    for (int i=1;i<tc.size();i+=2)
    {
        tc[i] = tsc.x*(tc[i] + tsh.x);
        tc[i-1] = tsc.y*(tc[i-1] + tsh.y);
    }
}
void TextureAtlas::process_tc(int num, glm::vec2 &tc)
{
    glm::vec2 tsc(1.0/gridWN,1.0/gridHN);
    glm::vec2 tsh(num % gridWN, num / gridWN);
    tc.x = tsc.x*(tc.x + tsh.x);
    tc.y = tsc.y*(tc.y + tsh.y);
}
bool TextureAtlas::bind()
{
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorTex.texture, 0);
    return true;
}
bool TextureAtlas::target(int num)
{
    glm::ivec2 tsh(num % gridWN, num / gridWN);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    //glViewport(tsh.x, tsh.y, gridWN, gridHN);
    glViewport(0, 0, width, height);

    return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
}
bool TextureAtlas::clear()
{
  target(0);
  glClearColor(clearColor.x,clearColor.y,clearColor.z,clearColor.a);
  glClear( GL_COLOR_BUFFER_BIT );
}
glm::mat4 TextureAtlas::tex_transform(int num)
{
  float a = num % gridWN;
  a /= gridWN;
  float b = num / gridWN;
  b /= gridHN;
  float h1 = 1.0f / gridWN;
  float h2 = 1.0f / gridHN;
  glm::mat4 sc_mat = glm::scale(glm::mat4(1.0f),glm::vec3(h1,h2,1));
  glm::mat4 tr_mat = glm::translate(glm::mat4(1.0f),glm::vec3(a,b,0));
  return tr_mat*sc_mat;
}
void TextureAtlas::gen_mipmaps()
{
    //do not work
    Model bm;
    std::vector<float> vertexes = {0,0,0,0,1,0,1,0,0,1,1,0};
    std::vector<float> tc = {0,0,0,0, 0,1,0,0, 1,0,0,0, 1,1,0,0};
    std::vector<GLuint> indices = {0,1,2,2,1,3};

    std::function<void(Model*)> _c_mip = [&](Model* h)
    {
        bm.positions = vertexes;
        bm.colors = tc;
        bm.indices = indices;
    };
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, width, height);
    glBindTexture(GL_TEXTURE_2D, colorTex.texture);
    int w = width;
    int h = height;
    int level = 1;
    while (w>1 && h>1)
    {   glBindTexture(GL_TEXTURE_2D, colorTex.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, level-1);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, level-1);
        Texture ctex(w,h,false);
        glBindTexture(GL_TEXTURE_2D, ctex.texture);
        glViewport(0, 0, w, h);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, ctex.texture, 0);
        glClearColor(clearColor.x,clearColor.y,clearColor.z,clearColor.a);
        glClear( GL_COLOR_BUFFER_BIT );

        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex",colorTex);
        bm.render(GL_TRIANGLES);

        glBindTexture(GL_TEXTURE_2D, ctex.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glDisable(GL_DEPTH_TEST);	
        glViewport(0, 0, w/2, h/2);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorTex.texture, level);
        glClearColor(clearColor.x,clearColor.y,clearColor.z,clearColor.a);
        glClear( GL_COLOR_BUFFER_BIT );
        bm.construct(_c_mip);
        mipMapRenderer.use();
        mipMapRenderer.texture("tex",ctex);
        mipMapRenderer.uniform("screen_size",glm::vec4(w,h,0,0));
        bm.render(GL_TRIANGLES);
        glEnable(GL_DEPTH_TEST);	
        w /= 2;
        h /= 2;
        level++;
    }
    glBindTexture(GL_TEXTURE_2D, colorTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, level-1);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorTex.texture, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}