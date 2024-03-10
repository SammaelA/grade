#include "texture_atlas.h"
#include <iostream>
#include "common_utils/matrix_transform.h"
#include "tinyEngine/model.h"
#include "tinyEngine/engine.h"
//long TextureAtlas::count = 0;
int atlases_count = 0;

TextureAtlas::TextureAtlas(): 
                              Countable(1),
                              colorTex(engine::textureManager->empty()),
                              normalTex(engine::textureManager->empty()),
                              mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
                              copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{

    fbo = create_framebuffer();
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}
TextureAtlas::TextureAtlas(const TextureAtlas &atlas):
Countable(1),
colorTex(atlas.colorTex),
normalTex(atlas.normalTex),
mipMapRenderer(atlas.mipMapRenderer),
copy(atlas.copy)
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
    fbo = create_framebuffer();
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(0,0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}
TextureAtlas::TextureAtlas(int w, int h, int l, int mip_levels, bool need_normals) :
                           Countable(1),
                           colorTex(engine::textureManager->create_texture_array(w, h, l, GL_RGBA8, mip_levels)),
                           normalTex(need_normals ? engine::textureManager->create_texture_array(w, h,l, GL_RGBA8, mip_levels) : 
                                                    Texture()),
                           mipMapRenderer({"mipmap_render.vs", "mipmap_render.fs"}, {"in_Position", "in_Tex"}),
                           copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"})
{
    debugl(10, "atlas created %dx%dx%d %f Mbytes\n",w,h,l,1e-6*2*w*h*l*4);

    width = w;
    height = h;
    layers = l;
    valid = true;
    fbo = create_framebuffer();
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(0,0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}
TextureAtlas::~TextureAtlas()
{
    //atlases_count--;
    debugl(10, "atlas deleted %dx%dx%d %f Mbytes\n",width,height,layers,1e-6*2*width*height*layers*4);
    delete_framebuffer(fbo);
}
void TextureAtlas::destroy()
{
    if (valid)
    {
        engine::textureManager->delete_tex(colorTex);
        engine::textureManager->delete_tex(normalTex);
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
float4 TextureAtlas::tc_transform(int num) const
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    float2 tsc(1.0 / gridWN, 1.0 / gridHN);
    float2 tsh(num % gridWN + l, num / gridWN);

    return float4(tsc.x,tsc.y,tsh.x,tsh.y);
}
void TextureAtlas::process_tc(int num, float3 &tc) const
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    float2 tsc(1.0 / gridWN, 1.0 / gridHN);
    float2 tsh(num % gridWN, num / gridWN);
    tc.x = tsc.x * (tc.x + tsh.x);
    tc.y = tsc.y * (tc.y + tsh.y);
    tc.z = l;
}
void TextureAtlas::pixel_offsets(int num, int3 &pix_tc) const
{
    float3 tc = float3(0,0,0);
    process_tc(num, tc);
    pix_tc = int3(width*tc.x, height*tc.y, tc.z);
}
void TextureAtlas::pixel_offsets(int num, uint4 &pix_tc) const
{
    float3 tc = float3(0,0,0);
    process_tc(num, tc);
    pix_tc = uint4(width*tc.x, height*tc.y, tc.z, 1);
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

    int2 tsh(num % gridWN, num / gridWN);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(l, type);
    glViewport(0, 0, width, height);

    return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
}
bool TextureAtlas::target_slice(int num, int type)
{     
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    int2 tsh(num % gridWN, num / gridWN);
    int2 sz(width/gridWN, height/gridHN);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    bind(l, type);
    glViewport(tsh.x*sz.x, tsh.y*sz.y, sz.x, sz.y);

    return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
}
bool TextureAtlas::clear()
{
    target(0,0);
    glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.w);
    glClear(GL_COLOR_BUFFER_BIT);

    target(0,1);
    glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.w);
    glClear(GL_COLOR_BUFFER_BIT);

    return true;
}
float4x4 TextureAtlas::tex_transform(int num) const
{
    int l = num / (gridHN*gridWN);
    num = num % (gridHN*gridWN);

    float a = num % gridWN;
    a /= gridWN;
    float b = num / gridWN;
    b /= gridHN;
    float h1 = 1.0f / gridWN;
    float h2 = 1.0f / gridHN;
    float4x4 sc_mat = LiteMath::scale(float4x4(), float3(h1, h2, 1));
    float4x4 tr_mat = LiteMath::translate(float4x4(), float3(a, b, 0));
    return tr_mat * sc_mat;
}
void TextureAtlas::gen_mipmaps(std::string mipmap_shader_name)
{
    debugl(10,"generate mipmaps\n");
    //logerr("mipmaps gen %d", tex(0).get_mip_levels());
    Shader copy({"copy_arr.vs", "copy_arr.fs"}, {"in_Position", "in_Tex"});
    Shader mipMapRenderer({"mipmap_render.vs", mipmap_shader_name}, {"in_Position", "in_Tex"});
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
            int mips = tex(k).get_mip_levels();

            for (int i=1;i<mips;i++)
            {
                glBindTexture(tex(k).type, tex(k).texture);
                glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BASE_LEVEL, i - 1);
                glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_LEVEL, i - 1);
                Texture ctex(engine::textureManager->create_texture(w,h));
                fbo1 = create_framebuffer();
                int prev_FBO = 0;
                glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
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
                mipMapRenderer.uniform("screen_size", float4(w, h, 0, 0));
                glDisable(GL_DEPTH_TEST);
                bm.render(GL_TRIANGLES);
                glEnable(GL_DEPTH_TEST);
                glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
                delete_framebuffer(fbo1);
                w /= 2;
                h /= 2;
                engine::textureManager->delete_tex(ctex);
            }
        }
        glBindTexture(tex(k).type, tex(k).texture);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_LEVEL, tex(k).get_mip_levels());
    }

    debugl(10,"generate mipmaps end\n");
}
TextureAtlas &TextureAtlas::operator=(const TextureAtlas &atlas)
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

    return *this;
}
Texture &TextureAtlas::tex(int type) 
{
    if (type == 0)
        return colorTex;
    else if (type == 1)
        return normalTex;
    
    return colorTex;
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
    Texture new_tex = engine::textureManager->create_texture_array(width, height, new_layers, GL_RGBA8, t.get_mip_levels());
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
    
    engine::textureManager->delete_tex(t);
    t = new_tex;
}

std::vector<int> TextureAtlas::get_all_valid_slices_ids()
{
    std::vector<int> ids;
    for (int i=0;i<curNum;i++)
    {
        ids.push_back(i);
    }
    return ids;
}

int TextureAtlasRawData::get_pixel_pos(int ww, int hh, int tex_id)
{
    int l = tex_id / (gridHN*gridWN);
    tex_id = tex_id % (gridHN*gridWN);
    int grid_h = tex_id / gridWN;
    int grid_w = tex_id % gridWN;

    return 4*w*h*l + 4*((hh + grid_h*slice_h)*w + (ww + grid_w*slice_w));
}
unsigned char TextureAtlasRawData::get_pixel_uc(int ww, int hh, Channel chan, int tex_id)
{

    int id =  get_pixel_pos(ww, hh, tex_id) + (int)chan;

    return raw_data[id];
}
unsigned char TextureAtlasRawData::get_pixel_uc_safe(int ww, int hh, Channel chan, int tex_id)
{
    if (!valid || !raw_data || ww < 0 || ww >= w || hh < 0 || hh >=h || tex_id < 0 || tex_id >=slices)
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
TextureAtlasRawData::TextureAtlasRawData(const TextureAtlas &atlas)
{
    if (!atlas.is_valid())
    {
        logerr("error invalid atlas for TextureAtlasRawData");
    }
    auto sizes = atlas.get_sizes();
    w = sizes.x;
    h = sizes.y;
    gridWN = sizes.z;
    gridHN = sizes.w;
    slice_h = h/gridHN;
    slice_w = w/gridWN;
    layers = atlas.layers_count();
    slices = layers * gridWN * gridHN;
    raw_data = new unsigned char[4*w*h*layers+1];

  /*int texDims[10];
    glBindTexture(GL_TEXTURE_2D_ARRAY, atlas.colorTex.texture);
    glGetTexLevelParameteriv(GL_TEXTURE_2D_ARRAY, 0, GL_TEXTURE_WIDTH, texDims);
    glGetTexLevelParameteriv(GL_TEXTURE_2D_ARRAY, 0, GL_TEXTURE_HEIGHT, texDims + 1);
    glGetTexLevelParameteriv(GL_TEXTURE_2D_ARRAY, 0, GL_TEXTURE_DEPTH, texDims + 2);
    logerr("tex sizes real %d %d %d but layers %d",texDims[0], texDims[1], texDims[2], layers);*/
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
    if (raw_data)
        delete[] raw_data;
    raw_data = nullptr;
    valid = false;
}
TextureAtlasRawData::~TextureAtlasRawData()
{
    clear();
}
int TextureAtlasRawData::get_slice_size(int tex_id)
{
    return 4*slice_h*slice_w;
}
void TextureAtlasRawData::get_slice(int tex_id, unsigned char *data, int *sw, int *sh)
{
    int ps = 0;
    for (int i=0;i<slice_h;i++)
    {
        memcpy(data + ps, raw_data + get_pixel_pos(0,i,tex_id), 4*slice_w);
        ps+=4*slice_w;
    }
    if (sw)
        *sw = slice_w;
    if (sh)
        *sh = slice_h;
}