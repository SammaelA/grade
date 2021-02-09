#pragma once
#include <map>
#include "tinyEngine/utility/texture.h"
#include <string>

class TextureManager
{
public:
    Texture get(std::string name);
    Texture get(int n);
    Texture get_arr(int n);
    Texture create_unnamed(int w, int h, bool shadow = false);
    Texture create_unnamed_array(int w, int h, bool shadow, int layers);
    Texture create_unnamed(SDL_Surface *s);
    Texture get_unnamed(GLuint n);
    Texture get_unnamed_arr(GLuint n);
    Texture load_unnamed(Texture &stub, unsigned char *data);
    Texture load_unnamed_arr(Texture &stub, unsigned char *data);
    Texture empty();
    TextureManager();
    TextureManager(std::string base_path);
    bool is_correct(Texture &t);
    void clear_unnamed();

private:
    std::map<std::string, Texture> textures;
    std::map<GLuint, Texture> unnamed_textures;
    std::map<GLuint, Texture> unnamed_array_textures;
};
extern TextureManager textureManager;