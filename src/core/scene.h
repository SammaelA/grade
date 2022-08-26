#pragma once
#include "grove.h"
#include "graphics_utils/terrain.h"
#include "grass.h"
#include "common_utils/bvh.h"

struct TreeTypeData;
struct Scene
{
    friend class boost::serialization::access;
    
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
        friend class boost::serialization::access;

        std::string name;
        ComplexModel model;
        std::vector<glm::mat4> instances;
        std::vector<AABB> bboxes;
        InstancedModel();

      private:
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & name;
          ar & model;
          ar & instances;
          ar & bboxes;
        }
    };
    Scene(){};
    void clear();
    ~Scene() {clear();}
    Scene& operator=(Scene &s) = delete;
    Scene& operator=(Scene &&s) = delete;
    
    Heightmap *heightmap = nullptr;
    std::vector<InstancedModel> instanced_models;

    GrovePacked grove;
    GrassPacked grass;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & heightmap;
      ar & instanced_models;
      ar & grove;
      ar & grass;
    }
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