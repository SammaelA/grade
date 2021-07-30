#include "texture_atlas.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "tinyEngine/utility/model.h"
#include "texture_manager.h"
//long TextureAtlas::count = 0;
int atlases_count = 0;

TextureAtlas::TextureAtlas(): colorTex(textureManager.empty()),
                              normalTex(textureManager.empty()),
                              mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
                              copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
TextureAtlas::TextureAtlas(const TextureAtlas &atlas):
mipMapRenderer(atlas.mipMapRenderer),
copy(atlas.copy),
colorTex(atlas.colorTex),
normalTex(atlas.normalTex)
{
    curNum = atlas.curNum;
    width = atlas.width;
    height = atlas.height;
    gridWN = atlas.gridWN;
    gridHN = atlas.gridHN;
    isGrid = atlas.isGrid;
    clearColor = atlas.clearColor;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(0,0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
TextureAtlas::TextureAtlas(int w, int h, int l) : 
                           colorTex(textureManager.create_unnamed_array(w, h, false, l)),
                           normalTex(textureManager.create_unnamed_array(w, h, false, l)),
                           mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
                           copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{
    logerr("atlas created %dx%d %f Mbytes",w,h,1e-6*2*w*h*4);
    //atlases_count++;
    width = w;
    height = h;
    layers = l;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(0,0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
TextureAtlas::~TextureAtlas()
{
    //atlases_count--;
    glDeleteFramebuffers(1, &fbo);
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
    if (curNum > layers * gridWN * gridHN)
        return -1;
    else
        return curNum - 1;
}
glm::vec4 TextureAtlas::tc_transform(int num)
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    glm::vec2 tsc(1.0 / gridWN, 1.0 / gridHN);
    glm::vec2 tsh(num % gridWN + l, num / gridWN);

    return glm::vec4(tsc.x,tsc.y,tsh.x,tsh.y);
}
void TextureAtlas::process_tc(int num, glm::vec3 &tc)
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    glm::vec2 tsc(1.0 / gridWN, 1.0 / gridHN);
    glm::vec2 tsh(num % gridWN, num / gridWN);
    tc.x = tsc.x * (tc.x + tsh.x);
    tc.y = tsc.y * (tc.y + tsh.y);
    tc.z = l;
}
bool TextureAtlas::bind(int layer, int type)
{
    if (type == 0)
        glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorTex.texture, 0, layer);
    else if (type == 1)
        glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, normalTex.texture, 0, layer);
    return true;
}
bool TextureAtlas::target(int num, int type)
{   
    
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    glm::ivec2 tsh(num % gridWN, num / gridWN);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(l, type);
    glViewport(0, 0, width, height);

    return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
}
bool TextureAtlas::clear()
{
    target(0,0);
    glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT);

    target(0,1);
    glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT);
}
glm::mat4 TextureAtlas::tex_transform(int num)
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    float a = num % gridWN;
    a /= gridWN;
    float b = num / gridWN;
    b /= gridHN;
    float h1 = 1.0f / gridWN;
    float h2 = 1.0f / gridHN;
    glm::mat4 sc_mat = glm::scale(glm::mat4(1.0f), glm::vec3(h1, h2, 1));
    glm::mat4 tr_mat = glm::translate(glm::mat4(1.0f), glm::vec3(a, b, 0));
    return tr_mat * sc_mat;
}
void TextureAtlas::gen_mipmaps()
{
    Shader copy({"copy_arr.vs", "copy_arr.fs"}, {"in_Position", "in_Tex"});
    Shader mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"});
    Model bm;
    std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
    std::vector<float> tc = {0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0};
    std::vector<GLuint> indices = {0, 1, 2, 2, 1, 3};

    std::function<void(Model *)> _c_mip = [&](Model *h) {
        bm.positions = vertexes;
        bm.colors = tc;
        bm.indices = indices;
    };

    for (int k=0;k<tex_count();k++)
    {
        GLuint fbo1;
        for (int j=0;j<layers;j++)
        {
            int w = get_sizes().x;
            int h =  get_sizes().y;
            int mips = 4;
            for (int i=1;i<mips;i++)
            {
                glBindTexture(tex(k).type, tex(k).texture);
                glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BASE_LEVEL, i - 1);
                glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_LEVEL, i - 1);
                Texture ctex(textureManager.create_unnamed(w,h));
                glGenFramebuffers(1, &fbo1);
                glBindFramebuffer(GL_FRAMEBUFFER, fbo1);
                glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, ctex.texture, 0);
                glViewport(0, 0, w, h);
                glClearColor(0,0,0,0);
                glClear(GL_COLOR_BUFFER_BIT);
                    bm.construct(_c_mip);
                    copy.use();
                    copy.texture("tex", tex(k));
                    copy.uniform("layer",(float)j);
                    bm.render(GL_TRIANGLES);
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                glBindTexture(ctex.type, ctex.texture);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
                glBindFramebuffer(GL_FRAMEBUFFER, fbo);
                glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, tex(k).texture, i, j);
                glViewport(0, 0, w / 2, h / 2);
                glClearColor(0,0,0,0);
                glClear(GL_COLOR_BUFFER_BIT);
                bm.construct(_c_mip);
                mipMapRenderer.use();
                mipMapRenderer.texture("tex", ctex);
                mipMapRenderer.uniform("screen_size", glm::vec4(w, h, 0, 0));
                glDisable(GL_DEPTH_TEST);
                bm.render(GL_TRIANGLES);
                glEnable(GL_DEPTH_TEST);
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                glDeleteFramebuffers(1, &fbo1);
                w /= 2;
                h /= 2;
            }
        }
    }
}
TextureAtlas &TextureAtlas::operator=(TextureAtlas &atlas)
{
    curNum = atlas.curNum;
    width = atlas.width;
    height = atlas.height;
    layers = atlas.layers;
    gridWN = atlas.gridWN;
    gridHN = atlas.gridHN;
    isGrid = atlas.isGrid;
    clearColor = atlas.clearColor;
    colorTex = atlas.colorTex;
    normalTex = atlas.normalTex;
    bind(0,0);
}
Texture &TextureAtlas::tex(int type)
{
    if (type == 0)
        return colorTex;
    else if (type == 1)
        return normalTex;
}