#pragma once 
#include "tinyEngine/resources.h"
#include <SDL2/SDL.h>
#include <functional>
#include "common_utils/utility.h"
class TextureManager;
class Texture
{
  friend class TextureManager;
  static constexpr GLuint INVALID_ID = ~0u;
public:
  GLuint texture = INVALID_ID;  //Texture int (OpenGL: everything is an int!)
  GLenum type = GL_TEXTURE_2D;  //Texture type (default is GL_TEXTURE_2D)
  GLenum format = GL_RGBA8;
  std::string origin = "";//empty for runtime generated textures, path to file for loaded textures
  int W = 0,H = 0,layers = 0;
  int tag = 0;//tag show to which part of program this texture belongs. Used by texture manager
  int mip_levels = 0;
  ~Texture();
  Texture();
  Texture(const Texture &t);
  Texture &operator=(const Texture &t);
  bool is_valid() const {return texture != INVALID_ID && W > 1 && H > 1 && layers > 0;}
  bool is_array() const {return is_valid() && layers>1;}
  int get_W() { return W; }
  int get_H() { return H; }
  int get_layers() { return layers; }
  int get_mip_levels() { return mip_levels; }
protected:
  Texture(GLuint _texture, GLenum _type, int _W, int _H, int _layers, int _tag, int _mip_levels, 
          GLenum _format, std::string _origin = "");
};