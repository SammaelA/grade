#pragma once 
#include "tinyEngine/texture.h"
#include "tinyEngine/postfx.h"

class ImageArithmetics
{
public:
 static void add(Texture &I, Texture &I_1, Texture &I_2, float a, float b); //I = a*I_1 + b*I_2
 static void mul(Texture &I, Texture &I_1, Texture &I_2, float a); //I = a*I_1*I_2
 static void div(Texture &I, Texture &I_1, Texture &I_2, float a); //I = a*I_1/I_2
};
