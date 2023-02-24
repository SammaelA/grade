#include "tinyEngine/engine.h"
#include "tinyEngine/image.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/model.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image_write.h"
#include "common_utils/blk.h"
#include <exception>

long tex_mem = 0;
int tex_count = 0;
long calls = 0;

TextureManager *engine::textureManager = nullptr;

Texture TextureManager::get(std::string name)
{
    auto t_it = textures.find(name);
    if (t_it == textures.end())
        return Texture();
    else if (t_it->second.texture == Texture::INVALID_ID)
    {
      textures.at(name) = load_tex(t_it->first, t_it->second.origin);
      return textures.at(name);
    }
    else
      return textures.at(name);
}
TextureManager::TextureManager()
{
    Texture t = empty();
    textures.emplace("empty",t);
}
TextureManager::~TextureManager()
{
  clear_unnamed();
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
    fbo1 = create_framebuffer();
    glBindFramebuffer(GL_FRAMEBUFFER, fbo1);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, t.texture, 0);
    for (int i=1;i<mips;i++)
    {
        glBindTexture(t.type, t.texture);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BASE_LEVEL, i - 1);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_LEVEL, i - 1);
        Texture ctex(engine::textureManager->create_texture(w,h));
        fbo = create_framebuffer();
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
        delete_framebuffer(fbo);
        w /= 2;
        h /= 2;
    }
    delete_framebuffer(fbo1);
}
TextureManager::TextureManager(std::string base_path, Block &textures_used)
{
    image::base_img_path = base_path;
    std::vector<std::string> names = {"texture not found"};
    std::vector<std::string> paths = {"wood3.jpg"};
    
    for (int i=0;i<textures_used.size();i++)
    {
        std::string name = textures_used.get_name(i);
        std::string path = textures_used.get_string(i,"");
        if (path != "")
        {
            names.push_back(name);
            paths.push_back(path);
        }
    }
    
    for (int i=0;i<paths.size();i++)
    {
      textures.emplace(names[i], Texture(Texture::INVALID_ID, GL_TEXTURE_2D, 1,1,1,0,1,GL_RGBA8,image::base_img_path + paths[i]));
    }
    debugl(10,"textures loaded %d\n",textures.size());
}
bool TextureManager::load_tex_to_catalog(std::string name, std::string path)
{
  if (get(name).texture != Texture::INVALID_ID)
    return true;
  
  Texture t = load_tex(name, path);
  if (t.texture != Texture::INVALID_ID)
  {
    textures.emplace(name, t);
    return true;
  }
  else
  {
    return false;
  }
}
Texture TextureManager::load_tex(std::string name, std::string path)
{
    try
    {
        auto ptr = image::load(path);
        if (!ptr)
            return Texture();
        Texture t = create_texture(ptr->w, ptr->h, GL_RGBA8, 9, ptr->pixels, GL_RGBA, GL_UNSIGNED_BYTE, path);
        mipmap(t, ptr->w, ptr->h, 9);
        SDL_FreeSurface(ptr);
        return t;
    }
    catch (const std::exception &e)
    {
        logerr("texture not found %s", path.c_str());
        return Texture();
    }
}

Texture TextureManager::load_unnamed_tex(std::string path, int mip_levels)
{
    try
    {
        auto ptr = image::load(path);
        if (!ptr)
            return empty();
        Texture t = create_texture(ptr->w, ptr->h, GL_RGBA8, mip_levels, ptr->pixels, GL_RGBA, GL_UNSIGNED_BYTE, path);
        mipmap(t, ptr->w, ptr->h, mip_levels);
        SDL_FreeSurface(ptr);
        unnamed_textures.emplace(t.texture, t);
        return t;
    }
    catch (const std::exception &e)
    {
        logerr("texture not found %s", path.c_str());
        return empty();
    }
}

Texture TextureManager::create_texture(int w, int h, GLenum format, int mip_levels, void *data,
                                       GLenum data_format, GLenum pixel_format, std::string origin_name)
{
  calls++;
  clean_unused();

  tex_count++;
  tex_mem += 4*w*h;
  GLuint texture;
  GLuint type = GL_TEXTURE_2D;
  glGenTextures(1, &texture);
  glBindTexture(type, texture);
  glTexStorage2D(texture, mip_levels, format, w, h);
  glTexImage2D(type, 0, format, w, h, 0, data_format, pixel_format, data);

  glTexParameteri(type, GL_TEXTURE_MIN_FILTER, mip_levels > 1 ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR);
  glTexParameteri(type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
  glTexParameteri(type, GL_TEXTURE_BASE_LEVEL, 0);
  if (mip_levels > 1)
    glTexParameteri(type, GL_GENERATE_MIPMAP, GL_TRUE);
  if (mip_levels > 1 && data)
    glGenerateMipmap(type);
  glBindTexture(type, 0);
  Texture t = Texture(texture, type, w, h, 1, current_textures_tag, mip_levels, format, origin_name);
  unnamed_textures.emplace(t.texture, t);
  return t;
}

Texture TextureManager::create_texture_cube(int w, int h, GLenum format, int mip_levels)
{
  calls++;
  clean_unused();

  tex_count++;
  tex_mem += 4*6*w*h;

  GLuint texture;
  GLuint type = GL_TEXTURE_CUBE_MAP;
  glGenTextures(1, &texture);
  glBindTexture(type, texture);
  glTexStorage2D(texture, mip_levels, format, w, h);

  glTexParameteri(type, GL_TEXTURE_MIN_FILTER, mip_levels > 1 ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR);
  glTexParameteri(type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
  glTexParameteri(type, GL_TEXTURE_BASE_LEVEL, 0);
  if (mip_levels > 1)
    glTexParameteri(type, GL_GENERATE_MIPMAP, GL_TRUE);
  glBindTexture(type, 0);
  Texture t = Texture(texture, type, w, h, 1, current_textures_tag, mip_levels, format, "");
  unnamed_cube_textures.emplace(t.texture, t);
  return t;
}

Texture TextureManager::create_texture_array(int w, int h, int layers, GLenum format, int mip_levels,
                                             void *data, GLenum data_format, GLenum pixel_format, std::string origin_name)
{
  calls++;
  clean_unused();

  tex_count++;
  tex_mem += 4*w*h*layers;

  GLuint texture;
  GLuint type = GL_TEXTURE_2D_ARRAY;
  glGenTextures(1, &texture);
  glBindTexture(type, texture);
  glTexImage3D(type, 0, GL_RGBA, w, h, layers, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
  int m_w = w/2;
  int m_h = h/2;
  for (int i=1;i<mip_levels;i++)
  {
    glTexImage3D(type, i, GL_RGBA, m_w, m_h, layers, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    m_w /= 2;
    m_h /= 2;
  }
  glTexParameteri(type, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
  glTexParameteri(type, GL_TEXTURE_BASE_LEVEL, 0);
  glTexParameteri(type, GL_TEXTURE_MAX_LEVEL, mip_levels - 1);
  if (mip_levels > 1)
    glTexParameteri(type, GL_GENERATE_MIPMAP, GL_TRUE);
  if (mip_levels > 1 && data)
    glGenerateMipmap(type);
  glBindTexture(type, 0);
  Texture t = Texture(texture, type, w, h, layers, current_textures_tag, mip_levels, format, origin_name);
  unnamed_array_textures.emplace(t.texture, t);
  return t;
}

Texture TextureManager::empty()
{
    return Texture(Texture::INVALID_ID, GL_TEXTURE_2D, 1, 1, 1, -1, 1, GL_RGBA8);
}
void TextureManager::clear_unnamed()
{
    for (auto const &p : unnamed_textures)
        glDeleteTextures(1, &(p.second.texture));
    
    for (auto const &p : unnamed_array_textures)
        glDeleteTextures(1, &(p.second.texture));
    
    for (auto const &p : unnamed_cube_textures)
        glDeleteTextures(1, &(p.second.texture));
}
void TextureManager::clean_unused()
{
  clear_unnamed_with_tag(-1);
}
void TextureManager::set_textures_tag(int tag)
{
    current_textures_tag = tag;
}
void TextureManager::clear_unnamed_with_tag(int tag)
{
    std::map<GLuint, Texture> unnamed_new;
    std::map<GLuint, Texture> unnamed_array_new;
    std::map<GLuint, Texture> unnamed_cube_new;

    for (auto const &p : unnamed_textures)
    {
        if ((p.second.tag != tag || tag < 0) && (ref_count.at(p.first) > 1))
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
        if ((p.second.tag != tag || tag < 0) && (ref_count.at(p.first) > 1))
            unnamed_array_new.emplace(p);
        else
        {
            tex_count--;
            tex_mem -= 4*p.second.W*p.second.H*p.second.layers;
            debugl(10,"unnamed texture array of size %dx%dx%d deleted\n",p.second.W, p.second.H,p.second.layers);
            glDeleteTextures(1, &(p.second.texture));
        }
    }

    for (auto const &p : unnamed_cube_textures)
    {
        if ((p.second.tag != tag || tag < 0) && (ref_count.at(p.first) > 1))
            unnamed_cube_new.emplace(p);
        else
        {
            tex_count--;
            tex_mem -= 4*6*p.second.W*p.second.H;
            debugl(10,"unnamed cubemap texture of size %dx%d deleted\n",p.second.W, p.second.H);
            glDeleteTextures(1, &(p.second.texture));
        }
    }

    unnamed_textures = unnamed_new;
    unnamed_array_textures = unnamed_array_new;
    unnamed_cube_textures = unnamed_cube_new;
    current_textures_tag = 0;
    debugl(10, "texture clearing completed. %d Mb allocated memory left\n", (int)(1e-6*tex_mem));
}
void TextureManager::save_png_raw(unsigned char *data, int w, int h, int channels, std::string name)
{
    save_png_raw_directly(data, w, h, channels, "saves/"+name+".png");
}
void TextureManager::save_png(Texture &t, std::string name)
{
    save_png_directly(t, "saves/"+name+".png");
}
void TextureManager::save_png_raw_directly(unsigned char *data, int w, int h, int channels, std::string name)
{
    if (data)
    {
        stbi_write_png(name.c_str(), w, h, channels, data, 0);
    }
    else
    {
        logerr("trying to save empty data");
    }
}
void TextureManager::save_png_directly(Texture &t, std::string name)
{
    unsigned char *data = nullptr;
    int w,h,layers = 1;
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
        layers = t.type == GL_TEXTURE_2D_ARRAY ? t.get_layers() : 1;
        data = safe_new<unsigned char>(4*w*h*layers, "save_png_data");

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
        save_png_raw_directly(data,w,h*layers,4,name);
        safe_delete<unsigned char>(data, "save_png_data"); 
    }
}

void TextureManager::save_bmp_raw(unsigned char *data, int w, int h, int channels, std::string name)
{
    save_bmp_raw_directly(data, w, h, channels, "saves/"+name+".bmp");
}
void TextureManager::save_bmp(Texture &t, std::string name)
{
    save_bmp_directly(t, "saves/"+name+".bmp");
}
void TextureManager::save_bmp_raw_directly(unsigned char *data, int w, int h, int channels, std::string name)
{
    if (data)
    {
        stbi_write_bmp(name.c_str(), w, h, channels, data);
    }
    else
    {
        logerr("trying to save empty data");
    }
}
void TextureManager::save_bmp_directly(Texture &t, std::string name)
{
    unsigned char *data = nullptr;
    int w,h,layers = 1;
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
        layers = t.type == GL_TEXTURE_2D_ARRAY ? t.get_layers() : 1;
        data = safe_new<unsigned char>(4*w*h*layers, "save_png_data");

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
        save_bmp_raw_directly(data,w,h*layers,4,name);
        safe_delete<unsigned char>(data, "save_png_data"); 
    }
}
void TextureManager::delete_tex(Texture &t)
{
  if (t.texture == Texture::INVALID_ID)
    return;
    if (t.type == GL_TEXTURE_2D)
    {
        auto it = unnamed_textures.find(t.texture);
        if (it == unnamed_textures.end())
        {
            logerr("trying to delete unregistered texture");
        }
        else
        {
            glDeleteTextures(1, &(it->second.texture));
            tex_count--;
            tex_mem -= 4*it->second.W*it->second.H*it->second.layers;
            unnamed_textures.erase(it);
        }
    }
    else if (t.type == GL_TEXTURE_2D_ARRAY)
    {
        auto it = unnamed_array_textures.find(t.texture);
        if (it == unnamed_array_textures.end())
        {
            logerr("trying to delete unregistered texture array");
        }
        else
        {
            glDeleteTextures(1, &(it->second.texture));
            tex_count--;
            tex_mem -= 4*it->second.W*it->second.H*it->second.layers;
            unnamed_array_textures.erase(it);
        }
    }
    else if (t.type == GL_TEXTURE_CUBE_MAP)
    {
        auto it = unnamed_cube_textures.find(t.texture);
        if (it == unnamed_cube_textures.end())
        {
            logerr("trying to delete unregisteredcubemap texture");
        }
        else
        {
            glDeleteTextures(1, &(it->second.texture));
            tex_count--;
            tex_mem -= 4*6*it->second.W*it->second.H*it->second.layers;
            unnamed_cube_textures.erase(it);
        }
    }
}