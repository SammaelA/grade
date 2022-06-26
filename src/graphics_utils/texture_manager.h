#pragma once
#include <map>
#include "../tinyEngine/texture.h"
#include <string>

struct Block;
class TextureManager
{
public:
    const static int baseMipLevelCount = 1;

    Texture get(std::string name);
    Texture create_texture(int w, int h, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount, void *data = nullptr,
                           GLenum data_format = GL_RGBA, GLenum pixel_format = GL_UNSIGNED_BYTE, std::string origin_name = "");
    Texture create_texture_cube(int w, int h, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount);
    Texture create_texture_array(int w, int h, int layers, GLenum format = GL_RGBA8, int mip_levels = baseMipLevelCount, 
                                 void *data = nullptr, GLenum data_format = GL_RGBA, GLenum pixel_format = GL_UNSIGNED_BYTE,
                                 std::string origin_name = "");
    Texture load_unnamed(Texture &stub, unsigned char *data);//obsolete, to be removed
    Texture load_unnamed_arr(Texture &stub, unsigned char *data);//obsolete, to be removed
    Texture empty();
    TextureManager();
    TextureManager(std::string base_path, Block &textures_used);
    bool load_tex(std::string name, std::string path);
    Texture load_unnamed_tex(std::string path);
    bool is_correct(Texture &t);
    void clear_unnamed();
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
    std::map<std::string, Texture> textures;
    std::map<GLuint, Texture> unnamed_textures;
    std::map<GLuint, Texture> unnamed_array_textures;
    std::map<GLuint, Texture> unnamed_cube_textures;

    int current_textures_tag = 0; 
};
extern TextureManager textureManager;