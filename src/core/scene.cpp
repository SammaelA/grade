#include "scene.h"
#include "graphics_utils/texture_manager.h"

Scene::~Scene() 
{
    if (heightmap) delete heightmap;
    for (auto &im : instanced_models)
        if (im.model)
            delete im.model;
}

Scene::InstancedModel::InstancedModel():
tex(textureManager.empty()),
model(nullptr)
{

}