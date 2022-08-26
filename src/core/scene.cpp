#include "scene.h"
#include "graphics_utils/texture_manager.h"

void Scene::clear() 
{
    if (heightmap) delete heightmap;
    heightmap = nullptr;
    for (auto &im : instanced_models)
      for (auto *m : im.model.models)
        if (m)
          delete m;
    instanced_models.clear();
}

Scene::InstancedModel::InstancedModel()
{

}