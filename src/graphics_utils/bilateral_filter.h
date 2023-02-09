#pragma once 
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"

class BilateralFilter
{
public:
 static Texture perform(Texture &tex, float sigma_d, float sigma_r);
};
