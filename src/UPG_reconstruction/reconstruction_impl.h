#pragma once
#include "reconstruction.h"
#include "tinyEngine/texture.h"
#include "generation_impl.h"
#include "tinyEngine/camera.h"

namespace upg
{
  //structure that contains all data for one given view, like mask and maybe some 
  //depth information
  struct ReferenceView
  {
    CameraSettings camera;
    bool fixed_camera = true;
    Texture mask;
    Texture resized_mask;
  };
}