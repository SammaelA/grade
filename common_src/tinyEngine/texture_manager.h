#pragma once
#include <map>
#include "texture.h"
#include <string>

struct Block;
class TextureManager
{
public:
  friend struct Texture;
    const static int baseMipLevelCount = 1;

    Texture get(std::string name);
    Texture create_texture(int w, int h, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount, void *data = nullptr,
                           GLenum data_format = GL_RGBA, GLenum pixel_format = GL_UNSIGNED_BYTE, std::string origin_name = "");
    Texture create_texture_cube(int w, int h, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount);
    Texture create_texture_array(int w, int h, int layers, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount, 
                                 void *data = nullptr, GLenum data_format = GL_RGBA, GLenum pixel_format = GL_UNSIGNED_BYTE,
                                 std::string origin_name = "");
    Texture empty();
    TextureManager();
    ~TextureManager();
    TextureManager(std::string base_path, Block &textures_used);
    bool load_tex_to_catalog(std::string name, std::string path);
    Texture load_unnamed_tex(std::string path, int mip_levels = 9);
    void set_textures_tag(int tag);
    void clear_unnamed_with_tag(int tag);
    void delete_tex(Texture &t);

    void save_png(Texture &t, std::string name);
    void save_png_raw(unsigned char *data, int w, int h, int channels, std::string name);
    void save_png_directly(Texture &t, std::string name);
    void save_png_raw_directly(unsigned char *data, int w, int h, int channels, std::string name);
    void save_bmp(Texture &t, std::string name);
    void save_bmp_raw(unsigned char *data, int w, int h, int channels, std::string name);
    void save_bmp_directly(Texture &t, std::string name);
    void save_bmp_raw_directly(unsigned char *data, int w, int h, int channels, std::string name);
private:
    void clear_unnamed();
    void clean_unused();
    Texture load_tex(std::string name, std::string path);

    std::map<GLuint, int> ref_count;
    std::map<std::string, Texture> textures;
    std::map<GLuint, Texture> unnamed_textures;
    std::map<GLuint, Texture> unnamed_array_textures;
    std::map<GLuint, Texture> unnamed_cube_textures;

    int current_textures_tag = 0; 
};