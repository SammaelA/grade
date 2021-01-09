#pragma once 
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <functional>
#include "../utility.h"

class TextureManager;
class Texture
{
  friend class TextureManager;
public:
  GLuint texture;               //Texture int (OpenGL: everything is an int!)
  GLenum type = GL_TEXTURE_2D;  //Texture type (default is this)
  void set_default_paramaters(Texture *tex);
  void bind();
  ~Texture();
  Texture(const Texture &t);
  Texture &operator=(const Texture &t);
protected:
  void clear();
  Texture();
  Texture(SDL_Surface* s);
  Texture(int W, int H, bool d = false);
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
