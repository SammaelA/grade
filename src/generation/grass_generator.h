#pragma once
#include "core/grass.h"
#include "common_utils/bbox.h"
#include "graphics_utils/terrain.h"
#include <map>
#include <vector>
#include <list>
#include "save_utils/serialization.h"

struct Cell;
struct Block;
struct LightVoxelsCube;
class GroveMask;
class GrassGenerator
{
public:
    friend class boost::serialization::access;

    void generate_grass_in_cell(Cell &cell, Field_2d *occlusion, GroveMask *mask, GroveMask *global_mask,
                                float density, std::vector<std::pair<int, float>> types);
    void remove_grass_in_cell(int cell_id);
    void pack_all_grass(GrassPacked &grass_packed, Heightmap &h);
private:
    struct GrassPatch
    {
        friend class boost::serialization::access;

        int grass_type;//position in used_grass_types vector
        Sphere2D sphere;
        int instances_max_cnt;
        std::list<Sphere2D> grass_instances;
      
    private:
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & grass_type;
          ar & sphere;
          ar & instances_max_cnt;
          ar & grass_instances;
        }
    };
    struct GrassInstance
    {
      friend class boost::serialization::access;

      glm::vec2 pos;
      float size;
      int cell_id;
      GrassInstance() {}
      GrassInstance(glm::vec2 _pos, float _size, int _cell_id)
      {
        pos = _pos;
        size = _size;
        cell_id = _cell_id;
      }

    private:
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar & pos;
        ar & size;
        ar & cell_id;
      }
    };
    std::vector<GrassPatch> grass_patches;
    std::vector<GrassType> used_grass_types;
    std::vector<float> grass_quantity;//for every type above
    std::vector<std::vector<GrassInstance>> grass_instances;
    void add_patch(std::vector<Cell> &cells, int cells_x, int cells_y, int cell_id, int type_id);

    //deprecated
    void prepare_grass_patches(std::vector<Cell> &cells, int cells_x, int cells_y);
    void set_grass_types(Block &grass_settings);
    void generate_grass_in_cell(Cell &cell, Field_2d *occlusion);

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & grass_patches;
      ar & used_grass_types;
      ar & grass_quantity;
      ar & grass_instances;
    }
};