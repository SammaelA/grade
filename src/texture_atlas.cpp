#include "texture_atlas.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "tinyEngine/utility/model.h"
#include "texture_manager.h"
//long TextureAtlas::count = 0;
int atlases_count = 0;

TextureAtlas::TextureAtlas(): 
                              Countable(1),
                              colorTex(textureManager.empty()),
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
Countable(1),
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
    occupied = atlas.occupied;
    valid = atlas.valid;
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
                           Countable(1),
                           colorTex(textureManager.create_unnamed_array(w, h, false, 4*l)),
                           normalTex(textureManager.create_unnamed_array(w, h, false, l)),
                           mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
                           copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{
    debugl(10, "atlas created %dx%d %f Mbytes\n",w,h,1e-6*2*w*h*4);

    width = w;
    height = h;
    layers = l;
    valid = true;
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
void TextureAtlas::destroy()
{
    if (valid)
    {
        textureManager.delete_tex(colorTex);
        textureManager.delete_tex(normalTex);
        valid = false;
        isGrid = false;
    }
}
void TextureAtlas::set_grid(int w, int h, bool _resizable)
{
    resizable = _resizable;
    curNum = 0;
    isGrid = true;
    gridWN = width / w;
    gridHN = height / h;
    occupied.resize(layers * gridWN * gridHN, false);
}
int TextureAtlas::add_tex()
{
    int sz = layers * gridWN * gridHN;
    if (curNum >= sz)
    {
        if (!resizable)
            return -1;
        else
        {
            increase_capacity();
        }
    }
    int pos = curNum;
    occupied.set(pos, true);
    curNum = sz;

    for (int i=pos+1;i<sz;i++)
    {
        if (!occupied.get_unsafe(i))
        {
            curNum = i;
            break;
        }
    }
    return pos;
}
void TextureAtlas::remove_tex(int pos)
{
    occupied.set(pos, false);
    curNum = MIN(curNum, pos);
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
    debugl(10,"generate mipmaps\n");
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
                textureManager.delete_tex(ctex);
            }
        }
    }
    debugl(10,"generate mipmaps end\n");
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
    valid = atlas.valid;
    bind(0,0);
}
Texture &TextureAtlas::tex(int type)
{
    if (type == 0)
        return colorTex;
    else if (type == 1)
        return normalTex;
}

int TextureAtlas::new_layers_count()
{
    return 1.25*layers + 1;
}
void TextureAtlas::increase_capacity()
{
    int new_layers = new_layers_count();
    int new_texs = (new_layers - layers)*gridHN*gridWN;

    for (int i=0;i<new_texs;i++)
    {
        occupied.push_back(false);
    }

    for (int k=0;k<tex_count();k++)
    {
        increase_capacity_tex(tex(k), new_layers);
    }
    debugl(10,"increase atlas size %d --> %d layers\n",layers, new_layers);
    layers = new_layers;
}
void TextureAtlas::increase_capacity_tex(Texture &t, int new_layers)
{
    Texture new_tex = textureManager.create_unnamed_array(width, height, false, new_layers);
    Shader copy({"copy_arr.vs", "copy_arr.fs"}, {"in_Position", "in_Tex"});
    Model bm;
    std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
    std::vector<float> tc = {0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0};
    std::vector<GLuint> indices = {0, 1, 2, 2, 1, 3};

    std::function<void(Model *)> _c_mip = [&](Model *h) {
        bm.positions = vertexes;
        bm.colors = tc;
        bm.indices = indices;
    };
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    for (int i=0;i<layers;i++)
    {
        glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, new_tex.texture, 0, i);
        glViewport(0, 0, width, height);
        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex", t);
        copy.uniform("layer",(float)i);
        bm.render(GL_TRIANGLES);
    }
    
    textureManager.delete_tex(t);
    t = new_tex;
}

unsigned char TextureAtlasRawData::get_pixel_uc(int ww, int hh, Channel chan, int tex_id)
{
    int grid_h = tex_id / gridWN;
    int grid_w = tex_id % gridWN;

    int id = 4*((hh + grid_h*slice_h)*w + (ww + grid_w*slice_w)) + (int)chan;

    return raw_data[id];
}
unsigned char TextureAtlasRawData::get_pixel_uc_safe(int ww, int hh, Channel chan, int tex_id)
{
    if (!valid || !raw_data || ww < 0 || w >= w || hh < 0 || hh >=h || tex_id < 0 || tex_id >=slices)
        return 0;
    else 
        return get_pixel_uc(ww,hh,chan,tex_id);
}
float TextureAtlasRawData::get_pixel(int w, int h, Channel chan, int tex_id)
{
    return (float)get_pixel_uc_safe(w, h, chan, tex_id)/255;
}
TextureAtlasRawData::TextureAtlasRawData()
{
    raw_data = nullptr;
    valid = false;
}
TextureAtlasRawData::TextureAtlasRawData(TextureAtlas &atlas)
{
    w = atlas.width;
    h = atlas.height;
    gridWN = atlas.gridWN;
    gridHN = atlas.gridHN;
    slice_h = h/gridHN;
    slice_w = w/gridWN;
    slices = atlas.layers * atlas.gridWN * atlas.gridHN;
    //raw_data = safe_new<unsigned char>(4*w*h*layers, "tex_atlas_raw_data");
    raw_data = new unsigned char[4*w*h*layers+1];
    glBindTexture(GL_TEXTURE_2D_ARRAY, atlas.colorTex.texture);

    glGetTexImage(GL_TEXTURE_2D_ARRAY,
                  0,
                  GL_RGBA,
                  GL_UNSIGNED_BYTE,
                  raw_data);
    glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

    valid = true;
}
void TextureAtlasRawData::clear()
{
    //safe_delete<unsigned char>(raw_data, "tex_atlas_raw_data");
    //if (raw_data)
    //    delete raw_data;
    raw_data = nullptr;
    valid = false;
}
TextureAtlasRawData::~TextureAtlasRawData()
{
    clear();
}