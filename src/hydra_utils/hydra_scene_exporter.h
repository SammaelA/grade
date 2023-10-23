#pragma once
#include <string>
struct Scene;
struct TreeTypeData;
struct Block;
namespace hydra 
{
  void get_default_settings(Block &b);
  bool export_scene(std::string directory, Scene &scene, Block &export_settings);
  bool export_internal(std::string directory, Scene &scene, Block &export_settings);
};