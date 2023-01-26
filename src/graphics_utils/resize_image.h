#pragma once
#include "tinyEngine/texture.h"
#include <vector>

class ImageResizer
{
public:
  ImageResizer();
  enum Type
  {
    STRETCH, //do not preserve x and y scale: 64x128 -> 64x64 will stretch original image in 2 times by x axis
    CENTERED //preserve x and y scale: 64x128 -> 64x64 will reduce the size of image by 2 and place in [16, 0]-[48, 64] box
  };
  static Texture resize(Texture tex, int new_w, int new_h, Type type, glm::vec4 base_color = glm::vec4(0,0,0,1));
};