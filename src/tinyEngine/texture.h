#pragma once 
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <functional>
#include "common_utils/utility.h"
#include "save_utils/saver.h"
class TextureManager;
class Texture
{
  friend class TextureManager;
  friend bool saver::save(FILE *f, Texture &t);
  friend bool saver::load(FILE *f, Texture &t);
  static constexpr GLuint INVALID_ID = ~0u;
public:
  GLuint texture = INVALID_ID;  //Texture int (OpenGL: everything is an int!)
  GLenum type = GL_TEXTURE_2D;  //Texture type (default is GL_TEXTURE_2D)
  GLenum format = GL_RGBA8;
  std::string origin = "";//empty for runtime generated textures, path to file for loaded textures
  int W = 0,H = 0,layers = 0;
  int tag = 0;//tag show to which part of program this texture belongs. Used by texture manager
  int mip_levels = 0;
  ~Texture()
  {

  }
  Texture()
  {

  }
  Texture(const Texture &t)
  {
    texture = t.texture;
    type = t.type;
    tag = t.tag;
    W = t.W;
    H = t.H;
    layers = t.layers;
    origin = t.origin;
    mip_levels = t.mip_levels;
    format = t.format;
  }
  Texture &operator=(const Texture &t)
  {
    texture = t.texture;
    type = t.type;
    tag = t.tag;
    W = t.W;
    H = t.H;
    layers = t.layers;
    origin = t.origin;
    mip_levels = t.mip_levels;
    format = t.format;

    return *this;
  }
  bool is_valid() {return W > 1 && H > 1 && layers > 0;}
  bool is_array() {return is_valid() && layers>0;}
  int get_W() { return W; }
  int get_H() { return H; }
  int get_layers() { return layers; }
  int get_mip_levels() { return mip_levels; }
protected:
  Texture(GLuint _texture, GLenum _type, int _W, int _H, int _layers, int _tag, int _mip_levels, 
          GLenum _format, std::string _origin = ""):
  texture(_texture),
  type(_type),
  W(_W),
  H(_H),
  layers(_layers),
  tag(_tag),
  mip_levels(_mip_levels),
  format(_format),
  origin(_origin)
  {
  };
};