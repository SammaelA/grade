#pragma once 
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"

class UnsharpMasking
{
public:
  static Texture perform(Texture &tex, float sigma = 1.5, float strength = 0.2);
};