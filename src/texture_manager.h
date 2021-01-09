#pragma once
#include <map>
#include "tinyEngine/utility/texture.h"
#include <string>

class TextureManager
{
public:
    Texture get(std::string name);
    Texture create_unnamed(int w, int h, bool shadow = false);
    Texture create_unnamed(SDL_Surface *s);
    Texture empty();
    TextureManager();
    TextureManager(std::string base_path);

private:
    std::map<std::string, Texture> textures;
    std::map<GLuint, Texture> unnamed_textures;
};
extern TextureManager textureManager;
