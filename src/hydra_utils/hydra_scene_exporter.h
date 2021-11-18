#pragma once
#include <string>
struct Scene;
struct TreeTypeData;
struct Block;
class HydraSceneExporter
{
public:
    bool export_scene(std::string directory, Scene &scene, Block &export_settings);
};