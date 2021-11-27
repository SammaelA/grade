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
public:
  GLuint texture = 0;               //Texture int (OpenGL: everything is an int!)
  GLenum type = GL_TEXTURE_2D;  //Texture type (default is this)
  std::string origin = "";//empty for runtime generated textures, path to file for loaded textures
  void set_default_paramaters(Texture *tex);
  void bind();
  ~Texture();
  Texture(const Texture &t);
  Texture &operator=(const Texture &t);
  bool is_valid() {return W > 1 && H > 1 && layers > 0;}
  bool is_array() {return is_valid() && layers>0;}
  int get_W() { return W; }
  int get_H() { return H; }
  int get_layers() { return layers; }
  int get_mip_levels() { return mip_levels; }
protected:
  int W,H,layers;
  int tag = 0;//tag show to which part of program this texture belongs. Used by texture manager
  unsigned id = 0;
  int mip_levels = 0;
  void clear();
  Texture();
  Texture(bool empty) {texture = ~0;}
  Texture(SDL_Surface* s);
  Texture(Texture &stub, unsigned char *data);
  Texture(int W, int H, bool d, int mip_levels);
  Texture(int W, int H, bool d, int layers, int mip_levels);
  void empty(int W, int H, bool set_default = true, GLenum F = GL_RGBA);
  void depth(int W, int H, bool set_default = true);
  void raw(SDL_Surface* s, bool set_default = true);
};
class Cubetexture: public Texture  //Cubetexture specialization.
{  
  friend class TextureManager;
public:                             //Same thing, 6 times
  Cubetexture(const Cubetexture &t);
  Cubetexture &operator=(const Cubetexture &t);
protected: 
  Cubetexture();
  Cubetexture(int W, int H, bool d = false);
  void empty(int W, int H, bool set_default = true, GLenum F = GL_RGBA);
  void depth(int W, int H, bool set_default = true);
};
