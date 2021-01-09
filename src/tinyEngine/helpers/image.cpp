#include "image.h"

namespace image {
  std::string base_img_path = "resources/textures/";
  SDL_Surface* load(std::string path){
    path = base_img_path + path;
    SDL_Surface* loaded = IMG_Load(path.c_str());
    if(loaded == nullptr)
    {
      printf( "Unable to load image %s! SDL_image Error: %s\n", path.c_str(), IMG_GetError() );
      return nullptr;
    }
    SDL_Surface* optimized = SDL_ConvertSurfaceFormat(loaded, SDL_PIXELFORMAT_RGBA32, 0);
    SDL_FreeSurface(loaded);
    return optimized;
  }

  void save(SDL_Surface* surface, std::string path){
    path = base_img_path + path;
    IMG_SavePNG(surface, path.c_str());
  }


}
