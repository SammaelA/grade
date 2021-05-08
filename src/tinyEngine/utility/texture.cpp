#include "texture.h"
        const int mipLevelCount = 5;
void Texture::set_default_paramaters(Texture *t)
{
        glTexParameteri(t->type, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(t->type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(t->type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(t->type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(t->type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
}
void Texture::bind()
{
        glBindTexture(type, texture);
}
Texture::Texture()
{
        glGenTextures(1, &texture);
        debugl(2,"texture created %d\n", texture);
} //Default constructor
Texture::Texture(SDL_Surface *s) : Texture()
{
        raw(s);
} //Load raw surface data into texture
Texture::Texture(int W, int H, bool d) : Texture()
{ //Create empty texture of defined size
        if (!d)
                empty(W, H);
        else
                depth(W, H);

        this->W = W;
        this->H = H;
        this->layers = 1;
}
Texture::Texture(Texture &stub, unsigned char *data) :
Texture()
{
        W = stub.W;
        H = stub.H;
        layers = stub.layers;
        type = stub.type;

        glBindTexture(type, texture);

        glTexParameteri(type, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
        glTexParameteri(type, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(type, GL_GENERATE_MIPMAP, GL_TRUE);
        if (type == GL_TEXTURE_2D)
        {
                glTexStorage2D(texture, mipLevelCount, GL_RGBA8, W, H);
                glTexImage2D(type, 0, GL_RGBA, W, H, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        }
        else 
        {
                glTexStorage3D(texture, mipLevelCount, GL_RGBA8, W, H, layers);
                glTexImage3D(type, 0, GL_RGBA, W, H,layers, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        }
        
        glGenerateMipmap(type);
}
Texture::~Texture()
{
        
}
//Texture Parameter Setter!
void Texture::empty(int W, int H, bool set_default, GLenum F)
{
        glBindTexture(type, texture);
        glTexImage2D(type, 0, F, W, H, 0, F, GL_UNSIGNED_BYTE, NULL);
        if (set_default)
                set_default_paramaters(this); //Call Parameter Setter
}

void Texture::depth(int W, int H, bool set_default)
{
        empty(W, H, set_default, GL_DEPTH_COMPONENT);
}

void Texture::raw(SDL_Surface *s, bool set_default)
{ //Generate a texture from raw surface data
        glBindTexture(type, texture);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        if (set_default)
                set_default_paramaters(this); //Call the parameter setting function!
        glTexImage2D(type, 0, GL_RGBA, s->w, s->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, s->pixels);
        glGenerateMipmap(GL_TEXTURE_2D);
        
        this->W = s->w;
        this->H = s->h;
        this->layers = 1;
}

  Texture::Texture(const Texture &t)
  {
          texture = t.texture;
          type = t.type;
          this->W = t.W;
          this->H = t.H;
          this->layers = t.layers;
  }
  Texture &Texture::operator=(const Texture &t)
  {
          texture = t.texture;
          type = t.type;
          this->W = t.W;
          this->H = t.H;
          this->layers = t.layers;
          return *this;
  }
  void Texture::clear()
  {
          glDeleteTextures(1, &texture);
  }
Texture::Texture(int W, int H, bool d, int layers)
{


        this->W = W;
        this->H = H;
        this->layers = layers;
        type = GL_TEXTURE_2D_ARRAY;
        glGenTextures(1, &texture);
        glBindTexture(type,texture);
        glTextureStorage3D(texture, mipLevelCount, GL_RGBA8, W, H, layers);
        glTexParameteri(type, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(type, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);
        glTexParameteri(type, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(type, GL_GENERATE_MIPMAP, GL_TRUE);
}
Cubetexture::Cubetexture() : Texture()
{
        type = GL_TEXTURE_CUBE_MAP;
}
Cubetexture::Cubetexture(int W, int H, bool d) : Cubetexture()
{ //Create empty texture of defined size
        if (!d)
                empty(W, H);
        else
                depth(W, H);
}

void Cubetexture::empty(int W, int H, bool set_default, GLenum F)
{
        glBindTexture(type, texture);
        for (unsigned int i = 0; i < 6; i++)
                glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, F, W, H, 0, F, GL_UNSIGNED_BYTE, NULL);
        if (set_default)
                set_default_paramaters(this);
}

void Cubetexture::depth(int W, int H, bool set_default)
{
        empty(W, H, set_default, GL_DEPTH_COMPONENT);
}