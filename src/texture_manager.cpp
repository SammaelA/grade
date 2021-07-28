#include "texture_manager.h"
#include "tinyEngine/helpers/image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "tinyEngine/helpers/stb_image_write.h"
#include <exception>

long tex_mem = 0;
int tex_count = 0;
Texture TextureManager::get_unnamed(GLuint n)
{
    auto t_it = unnamed_textures.find(n);
    if (t_it == unnamed_textures.end())
        return textures.at("texture not found");
    else
        return t_it->second;
}
Texture TextureManager::get_unnamed_arr(GLuint n)
{
    auto t_it = unnamed_array_textures.find(n);
    if (t_it == unnamed_array_textures.end())
        return textures.at("texture not found");
    else
        return t_it->second;
}
    Texture TextureManager::load_unnamed(Texture &stub, unsigned char *data)
    {
        Texture t = Texture(stub,data);
        unnamed_textures.emplace(t.texture,t);
        return t;
    }
    Texture TextureManager::load_unnamed_arr(Texture &stub, unsigned char *data)
    {
        Texture t = Texture(stub,data);
        unnamed_array_textures.emplace(t.texture,t);
        return t;
    }
bool TextureManager::is_correct(Texture &t)
{
    return t.texture != textures.at("texture not found").texture;
}
Texture TextureManager::get(std::string name)
{
    auto t_it = textures.find(name);
    if (t_it == textures.end())
        return textures.at("texture not found");
    else
        return textures.at(name);
}
Texture TextureManager::get(int n)
{
    n = n % textures.size();
    int i =0;
    for (auto it = textures.begin(); it != textures.end(); it++)
    {
        if (i == n)
            return it->second;
        i++;
    }
}
Texture TextureManager::get_arr(int n)
{
    n = n % unnamed_array_textures.size();
    int i =0;
    for (auto it = unnamed_array_textures.begin(); it != unnamed_array_textures.end(); it++)
    {
        if (i == n)
            return it->second;
        i++;
    }
}
TextureManager::TextureManager()
{
    Texture t(false);
    textures.emplace("empty",t);
}
void mipmap(Texture t, int w, int h, int mips = 9)
{
Shader copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"});
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
    GLuint fbo,fbo1;
    glGenFramebuffers(1, &fbo1);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo1);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, t.texture, 0);
    for (int i=1;i<mips;i++)
    {
        glBindTexture(t.type, t.texture);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BASE_LEVEL, i - 1);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_LEVEL, i - 1);
        Texture ctex(textureManager.create_unnamed(w,h));
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, ctex.texture, 0);
        glViewport(0, 0, w, h);
        glClearColor(1,0,0,1);
        glClear(GL_COLOR_BUFFER_BIT);
        
        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex", t);
        copy.uniform("layer",0);
        bm.render(GL_TRIANGLES);
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindTexture(ctex.type, ctex.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo1);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, t.texture, i);
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
        glDeleteFramebuffers(1, &fbo);
        w /= 2;
        h /= 2;
    }
    glDeleteFramebuffers(1, &fbo1);
}
TextureManager::TextureManager(std::string base_path)
{
    image::base_img_path = base_path;
    std::vector<std::string> names = {"reference_tree_test", "terrain","noise","colored_noise","leaf1","leaf2","leaf","wood","wood2","wood3",
                                      "grass","start_screen","texture not found"};
    std::vector<std::string> paths = {"reference_tree_test.png", "terrain.png","perlin.png","noise.png","leaf1.png","leaf2.png","leaf4.png","wood1.jpg","wood2.jpg","wood3.jpg",
                                      "grass.png","start_screen.png","texture_not_found.png"};
    for (int i=0;i<paths.size();i++)
    {
        try
        {
            auto ptr = image::load(image::base_img_path + paths[i]);
            if (!ptr)
                continue;
            Texture t(ptr);
            t.origin = paths[i];
            textures.emplace(names[i],t);
            mipmap(t,ptr->w,ptr->h,9);
            delete ptr;
        }
        catch(const std::exception& e)
        {
            logerr("texture not found %s",paths[i].c_str());
        }
    }
    debugl(10,"textures loaded %d\n",textures.size());
}

Texture TextureManager::create_unnamed(int w, int h, bool shadow)
{
    tex_mem += 4*w*h;
    tex_count++;
    debugl(10,"allocated %d Mb total for %d textures\n",(int)(1e-6*tex_mem), tex_count);
    Texture t(w,h,shadow);
    t.tag = current_textures_tag;
    unnamed_textures.emplace(t.texture,t);
    return t;
}
Texture TextureManager::create_unnamed_array(int w, int h, bool shadow, int layers)
{
    tex_mem += 4*w*h*layers;
    tex_count++;
    debugl(10,"allocated %d Mb total for %d textures\n",(int)(1e-6*tex_mem), tex_count);
    Texture t(w,h,shadow,layers);
    t.tag = current_textures_tag;
    unnamed_array_textures.emplace(t.texture,t);
    return t;
}
Texture TextureManager::create_unnamed(SDL_Surface *s)
{
    tex_mem += 4*s->w*s->h;
    tex_count++;
    debugl(10,"allocated %d Mb total for %d textures\n",(int)(1e-6*tex_mem), tex_count);
    Texture t(s);
    t.tag = current_textures_tag;
    unnamed_textures.emplace(t.texture,t);
    return t;
}
Texture TextureManager::empty()
{
    return get("empty");
}
void TextureManager::clear_unnamed()
{
    for (auto const &p : unnamed_textures)
        glDeleteTextures(1, &(p.second.texture));
    
    for (auto const &p : unnamed_array_textures)
        glDeleteTextures(1, &(p.second.texture));
}
void TextureManager::set_textures_tag(int tag)
{
    current_textures_tag = tag;
}
void TextureManager::clear_unnamed_with_tag(int tag)
{
    std::map<GLuint, Texture> unnamed_new;
    std::map<GLuint, Texture> unnamed_array_new;

    for (auto const &p : unnamed_textures)
    {
        if (p.second.tag != tag)
            unnamed_new.emplace(p);
        else
        {
            tex_count--;
            tex_mem -= 4*p.second.W*p.second.H;
            debugl(10,"unnamed texture of size %dx%d deleted\n",p.second.W, p.second.H);
            glDeleteTextures(1, &(p.second.texture));
        }
    }
    
    for (auto const &p : unnamed_array_textures)
    {
        if (p.second.tag != tag)
            unnamed_array_new.emplace(p);
        else
        {
            tex_count--;
            tex_mem -= 4*p.second.W*p.second.H*p.second.layers;
            debugl(10,"unnamed texture array of size %dx%dx%d deleted\n",p.second.W, p.second.H,p.second.layers);
            glDeleteTextures(1, &(p.second.texture));
        }
    }

    unnamed_textures = unnamed_new;
    unnamed_array_textures = unnamed_array_new;
    current_textures_tag = 0;
    debugl(10, "texture clearing completed. %d Mb allocated memory left\n", (int)(1e-6*tex_mem));
}

void TextureManager::save_bmp(Texture &t, std::string name)
{
    unsigned char *data = nullptr;
    int w,h;
    if (!t.is_valid())
    {
        logerr("trying to save invalid texture");
        data = nullptr;
        w = 0;
        h = 0;
    }
    else if (t.type == GL_TEXTURE_2D || t.type == GL_TEXTURE_2D_ARRAY)
    {
        w = t.get_W();
        h = t.get_H();
        data = safe_new<unsigned char>(4*w*h, "save_png_data");

        glBindTexture(t.type, t.texture);

        glGetTexImage(t.type,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    data);
        glBindTexture(t.type, 0);
    }
    else
    {
        logerr("invalid texture format");
        data = nullptr;
        w = 0;
        h = 0;
    }
    if (data)
    {
        std::string path = "saves/"+name+".bmp";
        stbi_write_bmp(path.c_str(), w, h, 4, data);
        safe_delete<unsigned char>(data, "save_png_data"); 
    }
}