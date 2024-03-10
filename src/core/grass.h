#pragma once
#include <vector>
#include "common_utils/LiteMath_ext.h"
#include "graphics_utils/texture_atlas.h"
#include "common_utils/blk.h"
#include "save_utils/serialization.h"

struct GrassType
{   
    friend class boost::serialization::access;
    
    int id = -1;
    Texture texture;
    float patch_size = 25;
    float patch_size_std_dev = 10;
    float patch_density = 1;
    float patch_density_std_dev = 0.5;
    float plant_size = 3;
    float plant_size_std_dev = 1.5;
    float light_sensivity = 1;
    float push = 1;
    std::string model_name = "simple_grass";
    void load_from_blk(Block &block);
    GrassType();

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & id;
      ar & texture;
      ar & patch_size;
      ar & patch_size_std_dev;
      ar & patch_density;
      ar & patch_density_std_dev;
      ar & plant_size;
      ar & plant_size_std_dev;
      ar & light_sensivity;
      ar & push;
      ar & model_name;
    }
};

struct GrassInstanceData
{   
    friend class boost::serialization::access;

    glm::vec3 pos;
    float size;
    float rot_y;
    GrassInstanceData() {};
    GrassInstanceData(glm::vec3 _pos, float _size, float _rot_y)
    {
        pos = _pos;
        size = _size;
        rot_y = _rot_y;
    }

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & pos;
      ar & size;
      ar & rot_y;
    }
};

struct GrassPacked
{   
    friend class boost::serialization::access;

    std::vector<GrassType> used_grass_types;
    TextureAtlas grass_textures;
    std::vector<std::pair<int, std::vector<GrassInstanceData>>> grass_instances;
    void clear()
    {
      used_grass_types.clear(); 
      grass_instances.clear();
    }
    ~GrassPacked() {clear();}
    //each member is a type: int is texture id in atlas, second is list of instances for this type

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & used_grass_types;
      ar & grass_textures;
      ar & grass_instances;
    }
};