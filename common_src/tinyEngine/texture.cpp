#include "texture.h"
#include "engine.h"

  Texture::~Texture()
  {
    if (texture != INVALID_ID)
    {
      auto it = engine::textureManager->ref_count.find(texture);
      if (it != engine::textureManager->ref_count.end())
        it->second--;
    }
  }
  Texture::Texture()
  {
    texture = INVALID_ID;
    type = GL_TEXTURE_2D;
  }
  Texture::Texture(const Texture &t)
  {
    if (t.texture != INVALID_ID)
    {
      engine::textureManager->ref_count.at(t.texture)++;
    }
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
  Texture &Texture::operator=(const Texture &t)
  {
    if (texture != INVALID_ID)
    {
      engine::textureManager->ref_count.at(texture)--;
    }
    if (t.texture != INVALID_ID)
    {
      engine::textureManager->ref_count.at(t.texture)++;
    }
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
  Texture::Texture(GLuint _texture, GLenum _type, int _W, int _H, int _layers, int _tag, int _mip_levels, 
          GLenum _format, std::string _origin):
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
    if (texture != INVALID_ID)
    {
      engine::textureManager->ref_count.emplace(texture, 0);//this is an initial texture, that is stored in texture manager,
                                                   //it is not counted in ref_count, as if we don't have references 
                                                   //to this texture outside of manager -> we don't need it
      engine::textureManager->ref_count.at(texture)++;
    }
  };