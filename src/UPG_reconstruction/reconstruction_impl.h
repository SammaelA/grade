#pragma once
#include "reconstruction.h"
#include "tinyEngine/texture.h"

namespace upg
{
  //structure that contains all data for one given view, like mask and maybe some 
  //depth information
  struct ReferenceView
  {
    Texture mask;
  };
}