#include "texture_manager.h"
#include "tinyEngine/helpers/image.h"
Texture &TextureManager::get(std::string name)
{
    auto t_it = textures.find(name);
    if (t_it == textures.end())
        return textures.at("texture not found");
    else
        return textures.at(name);
}
TextureManager::TextureManager(std::string base_path)
{
    image::base_img_path = base_path;

    textures.emplace("leaf", image::load("leaf_2.png"));
    textures.emplace("wood", image::load("bark-1.jpg"));
    textures.emplace("texture not found", image::load("leaf_2.png"));
}