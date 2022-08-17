#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
#include "grass.h"
#include "common_utils/bvh.h"

struct TreeTypeData;
struct Scene
{
    enum ObjCategories
    {
        UNKNOWN,
        ERROR,
        SIMPLE_OBJECT,
        DEBUG_MODEL,
        TREE
    };
    struct InstancedModel
    {
        std::string name;
        Texture tex;
        ComplexModel model;
        std::vector<glm::mat4> instances;
        std::vector<AABB> bboxes;
        InstancedModel();
    };
    Scene(){};
    ~Scene();
    Scene& operator=(Scene &&s)
    {
      if (heightmap)
        delete heightmap;
      heightmap = s.heightmap;
      s.heightmap = nullptr;

      for (auto &im : instanced_models)
        for (auto *m : im.model.models)
          if (m)
            delete m;
      instanced_models = std::move(instanced_models);

      grove = s.grove;
      grass = s.grass;

      return *this;
    }
    
    Heightmap *heightmap = nullptr;
    std::vector<InstancedModel> instanced_models;

    GrovePacked grove;
    GrassPacked grass;
};

  enum RenderPixelTypes
  {
    // should match with deffered_light.fs
    PIXEL_TYPE_NONE = 0,
    PIXEL_TYPE_TERRAIN = 1,
    PIXEL_TYPE_TREES = 2,
    PIXEL_TYPE_GRASS = 3,
    PIXEL_TYPE_MODELS = 4,
    PIXEL_TYPE_DEBUG_LIGHT = 5,
    PIXEL_TYPE_DEBUG_NO_LIGHT = 6
  };