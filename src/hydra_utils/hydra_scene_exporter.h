#pragma once
#include <string>
struct Scene;
struct TreeTypeData;
struct Block;
class HydraSceneExporter
{
public:
    bool export_scene(std::string directory, Scene &scene, Block &export_settings);
private:
    bool export_internal1(std::string directory, Scene &scene, Block &export_settings);
    bool export_internal2(std::string directory, Scene &scene, Block &export_settings);
};