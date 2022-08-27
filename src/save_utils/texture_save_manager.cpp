#include "texture_save_manager.h"
#include "tinyEngine/engine.h"

TextureSaveManager *texSaveManager = nullptr;

void TextureSaveManager::start_save(std::string _directory, bool _make_copies)
{
  directory = _directory;
  make_copies = _make_copies;
}
std::string TextureSaveManager::register_texture_to_save(const Texture &t)
{
  if (!t.is_valid())
    return "";
  else if (!make_copies && (t.origin != ""))
  {
    // this texture already exists on hard drive, no need to save
    return t.origin;
  }
  else
  {
    std::string type_s = "tex2D";
    if (t.type == GL_TEXTURE_2D)
      type_s = "tex2D";
    else if (t.type == GL_TEXTURE_2D_ARRAY)
      type_s = "tex3D";
    else if (t.type == GL_TEXTURE_CUBE_MAP)
      type_s = "texCube";
    else
      type_s = "texSpecial";
    std::string path = directory + "/" + type_s + "_" + std::to_string(t.texture) + ".tex";
    if (textures_to_save.find(path) == textures_to_save.end())
    {
      textures_to_save.emplace(path, t);
    }
    return path;
  }
}
void TextureSaveManager::finish_save()
{
  for (auto &p : textures_to_save)
  {
    Texture &t = p.second;

    unsigned char *pixels = new unsigned char[4 * t.W * t.H * t.layers];
    glBindTexture(t.type, t.texture);
    glGetTexImage(t.type, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    glBindTexture(t.type, 0);
    std::string name = p.first;
    FILE *tex = fopen(name.c_str(), "wb");
    fwrite(pixels, sizeof(unsigned char), 4 * t.W * t.H * t.layers, tex);
    fclose(tex);
    engine::textureManager->save_png_directly(t, name+".png");
    delete[] pixels;
  }
  textures_to_save.clear();
}
bool TextureSaveManager::load_texture(std::string token, Texture &t)
{
  if (token == "")
  {
    t = engine::textureManager->empty();
    return false;
  }
  if (token.substr(token.size() - 4, 4) == std::string(".tex"))
  {
    // raw texture
    unsigned char *pixels = new unsigned char[4 * t.W * t.H * t.layers];
    FILE *tex = fopen(token.c_str(), "rb");
    int c = fread(pixels, sizeof(unsigned char), 4 * t.W * t.H * t.layers, tex);

    if (t.type == GL_TEXTURE_2D)
      t = engine::textureManager->create_texture(t.W, t.H, t.format, t.mip_levels, (void *)pixels);
    else if (t.type == GL_TEXTURE_2D_ARRAY)
    {
      t = engine::textureManager->create_texture_array(t.W, t.H, t.layers, t.format, t.mip_levels, (void *)pixels);
    }
    else
    {
      t = engine::textureManager->empty();
    }
    fclose(tex);
    delete[] pixels;
    engine::textureManager->save_png_directly(t, token+"_loaded.png");
    //logerr("loaded texture %d", t.texture);
    return true;
  }
  else
  {
    t = engine::textureManager->load_unnamed_tex(token);
    return true;
  }
}