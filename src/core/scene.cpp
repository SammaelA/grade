#include "scene.h"
#include "graphics_utils/texture_manager.h"

Scene::~Scene() 
{
    if (heightmap) delete heightmap;
    for (auto &im : instanced_models)
      for (auto *m : im.model.models)
        if (m)
          delete m;
}

Scene::InstancedModel::InstancedModel():
tex(textureManager.empty())
{

}