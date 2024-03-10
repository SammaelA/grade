#pragma once
#include "core/scene.h"
#include "biome.h"
#include "grove_generation_utils.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "grove_generator.h"
#include "grove_packer.h"
#include "grass_generator.h"
#include <mutex>

struct SceneGenerationContext;
struct Cell;

namespace scene_gen
{
  struct Patch
  {
    friend class boost::serialization::access;

    Sphere2D border;
    float density;
    std::vector<std::pair<int, float>> types;
    int patch_type_id;
    int biome_id;
    Patch(){};
    Patch(float2 pos, Biome::PatchDesc &patchDesc, int biome_id, int patch_type_id);

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & border;
      ar & density;
      ar & types;
      ar & patch_type_id;
      ar & biome_id;
    }
  };

  void generate_plants_cells(SceneGenerationContext &ctx, std::vector<int> cell_ids);
  void generate_grass_cells(SceneGenerationContext &ctx, std::vector<int> cell_ids);
  void remove_trees_from_scene(SceneGenerationContext &ctx, std::vector<int> &ids);
  void prepare_tree_prototypes(SceneGenerationContext &ctx);
};

struct SceneGenerationContext
{
  friend class boost::serialization::access;

  SceneGenerationContext() : objects_bvh(false){};
  ~SceneGenerationContext() {if (inited) clear();}
  SceneGenerationContext &operator=(SceneGenerationContext &&ctx) = delete;
  Scene scene;
  Block settings;
  BVH objects_bvh;
  BiomeMap biome_map;
  GroveMask global_mask;
  GrovePacker packer;
  GrassGenerator grass_generator;
  std::vector<Cell> cells;
  std::vector<scene_gen::Patch> trees_patches;
  std::vector<scene_gen::Patch> grass_patches;
  int cells_x, cells_y;
  float hmap_pixel_size, biome_map_pixel_size;
  float2 heightmap_size, grass_field_size, full_size, cell_size, center, start_pos;
  float3 center3;

  bool inited = false;

  void clear()
  {
    scene.clear();
    settings = Block();
    objects_bvh.clear();
    biome_map = BiomeMap();
    global_mask = GroveMask();
    grass_generator = GrassGenerator();
    packer.clear();
    cells.clear();
    trees_patches.clear();
    grass_patches.clear();
  }

private:
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    {
      int b_id = AbstractTreeGenerator::branch_next_id;
      int t_id = AbstractTreeGenerator::tree_next_id;

      ar & b_id;
      ar & t_id;

      AbstractTreeGenerator::branch_next_id.store(b_id);
      AbstractTreeGenerator::tree_next_id.store(t_id);
    }
    

    ar & scene;
    if (cur_ser_helper)
      cur_ser_helper->instanced_branches.set_container(&(scene.grove.instancedBranches));
    ar & settings;
    ar & objects_bvh;
    ar & biome_map;
    ar & global_mask;
    //if (Archive::is_loading::value)
    //  packer.clear();
    ar & packer;
    ar & grass_generator;
    ar & cells;
    ar & trees_patches;
    ar & grass_patches;

    ar & cells_x;
    ar & cells_y;
    ar & hmap_pixel_size;
    ar & biome_map_pixel_size;
    ar & heightmap_size;
    ar & grass_field_size;
    ar & full_size;
    ar & cell_size;
    ar & center;
    ar & start_pos;
    ar & center3;
    ar & inited;
  }
};
struct Cell
{
  friend class boost::serialization::access;

  enum CellStatus
  {
    EMPTY,//nothing to do with it or task is not set
    WAITING,//task is set but not performed. voxels_small == nullptr, planar_occlusion == nullptr
    BORDER,//plants are generated. voxels_small!= nullptr, planar_occlusion == nullptr
    FINISHED_PLANTS, //plants are generated. voxels_small == nullptr, planar_occlusion != nullptr
    FINISHED_ALL //plants and grass are generated. voxels_small == nullptr, planar_occlusion == nullptr
  };
  std::mutex cell_lock;
  std::vector<GrovePrototype> prototypes;

  LightVoxelsCube *voxels_small = nullptr;//used for plants gen is dependant cells
  Field_2d *planar_occlusion = nullptr;//used for grass gen in this cell
  int id = -1;
  AABB2D bbox;
  std::atomic<CellStatus> status;
  std::vector<int> depends;
  std::vector<int> grass_patches;//set and used by grass generator
  std::vector<int> trees_patches;
  std::vector<std::pair<int,int>> biome_stat;//<biome_id, number_of_pixels> 
  AABB influence_bbox;
  explicit Cell(CellStatus _status = CellStatus::EMPTY) {status = _status;}
  Cell(const Cell &cell)
  {
    id = cell.id;
    prototypes = cell.prototypes;
    bbox = cell.bbox;
    status.store(cell.status);
    influence_bbox = cell.influence_bbox;
    depends = cell.depends;

    if (voxels_small || planar_occlusion || !grass_patches.empty())
    {
      logerr("only copying of empty cells is allowed");
    }
  };
  ~Cell()
  {
    if (voxels_small)
      delete voxels_small;
    if (planar_occlusion)
      delete planar_occlusion; 
    for (auto &prototype : prototypes)
    {
      if (prototype.biome_mask)
      {
        delete prototype.biome_mask;
        prototype.biome_mask = nullptr;
      }
    }
  }
private:
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & prototypes;
    ar & planar_occlusion;
    ar & id;
    ar & bbox;
    status = EMPTY;
    ar & depends;
    ar & grass_patches;
    ar & trees_patches;
    ar & biome_stat;
    ar & influence_bbox;
  }
};