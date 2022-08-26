#pragma once
#include "tinyEngine/texture.h"
#include <map>
class TextureSaveManager
{
public:
  void start_save(std::string directory, bool make_copies = false);
  std::string register_texture_to_save(const Texture &t);//returns token that can be saved and then 
                                                   //manager can restore this texture by token
  void finish_save();
  bool load_texture(std::string token, Texture &t);
private:
  std::map<std::string, Texture> textures_to_save;
  std::string directory = "saves";
  bool make_copies = false;
};

extern TextureSaveManager *texSaveManager;