#pragma once 
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <functional>
#include "../utility.h"
class Texture{
public:
  GLuint texture;               //Texture int (OpenGL: everything is an int!)
  GLenum type = GL_TEXTURE_2D;  //Texture type (default is this)
  void set_default_paramaters(Texture *tex);
  void bind() {glBindTexture( type, texture );}
  Texture(){ glGenTextures( 1, &texture ); debug("texture created %d\n",texture);};         //Default constructor
  Texture(SDL_Surface* s):Texture(){ raw(s); };     //Load raw surface data into texture
  Texture(int W, int H, bool d = false):Texture(){  //Create empty texture of defined size 
    if(!d) empty(W, H);
    else   depth(W, H);
  };

  ~Texture(){ glDeleteTextures(1, &texture); }
      //Texture Parameter Setter!
  void empty(int W, int H, bool set_default = true, GLenum F = GL_RGBA){
    glBindTexture( type, texture );
    glTexImage2D(type, 0, F, W, H, 0, F, GL_UNSIGNED_BYTE, NULL); 
    if (set_default) 
      set_default_paramaters(this);           //Call Parameter Setter
  }

  void depth(int W, int H, bool set_default = true){
    empty(W, H, set_default, GL_DEPTH_COMPONENT);
  }

  void raw(SDL_Surface* s, bool set_default = true){  //Generate a texture from raw surface data
    glBindTexture( type, texture );
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //glTexImage2D(type, 0, GL_RGBA, s->w, s->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, s->pixels);
    if (set_default) 
      set_default_paramaters(this);   //Call the parameter setting function!
    glTexImage2D(type, 0, GL_RGBA, s->w, s->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, s->pixels);
    glGenerateMipmap(GL_TEXTURE_2D);
  }
};

//Default parameter setting function!
//Note that you can pass your own to the texture building functions above!z


class Cubetexture: public Texture{  //Cubetexture specialization.
public:                             //Same thing, 6 times
  Cubetexture():Texture(){
    type = GL_TEXTURE_CUBE_MAP;
  };

  Cubetexture(int W, int H, bool d = false):Cubetexture(){  //Create empty texture of defined size
    if(!d) empty(W, H);
    else   depth(W, H);
  };

  void empty(int W, int H, bool set_default = true, GLenum F = GL_RGBA){
    glBindTexture(type, texture);
    for(unsigned int i = 0; i < 6; i++)
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, F, W, H, 0, F, GL_UNSIGNED_BYTE, NULL);
    if (set_default) 
      set_default_paramaters(this);  
  }

  void depth(int W, int H, bool set_default = true){
    empty(W, H, set_default, GL_DEPTH_COMPONENT);
  };
};
