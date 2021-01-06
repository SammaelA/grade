#pragma once
#include <map>
#include "tinyEngine/utility/texture.h"
#include <string>

class TextureManager
{
public:
    Texture &get(std::string name);
    TextureManager(){};
    TextureManager(std::string base_path);

private:
    std::map<std::string, Texture> textures;
};
extern TextureManager textureManager;
