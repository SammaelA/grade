#include "texture_manager.h"
#include "tinyEngine/helpers/image.h"
#include <exception>

Texture TextureManager::get(std::string name)
{
    auto t_it = textures.find(name);
    if (t_it == textures.end())
        return textures.at("texture not found");
    else
        return textures.at(name);
}
TextureManager::TextureManager()
{
    Texture t;
    textures.emplace("empty",t);
}
TextureManager::TextureManager(std::string base_path)
{
    image::base_img_path = base_path;
    std::vector<std::string> names = {"leaf","wood","texture not found"};
    std::vector<std::string> paths = {"leaf_2.png","bark-1.jpg","bark-1.jpg"};
    for (int i=0;i<3;i++)
    {
        try
        {
            auto ptr = image::load(image::base_img_path + paths[i]);
            if (!ptr)
                continue;
            Texture t(ptr);
            textures.emplace(names[i],t);
        }
        catch(const std::exception& e)
        {
            logerr("texture not found %s",paths[i].c_str());
        }
    }
    debug("textures loaded %d\n",textures.size());
}
Texture TextureManager::create_unnamed(int w, int h, bool shadow)
{
    Texture t(w,h,shadow);
    unnamed_textures.emplace(t.texture,t);
    return t;
}
Texture TextureManager::create_unnamed(SDL_Surface *s)
{
    Texture t(s);
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
}