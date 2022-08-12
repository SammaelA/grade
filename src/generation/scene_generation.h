#pragma once
#include "core/scene.h"
#include "biome.h"
#include "grove_generation_utils.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "grove_generator.h"
#include "grove_packer.h"
#include <mutex>

struct SceneGenerationContext;
struct Cell;

namespace scene_gen
{
  struct Patch
  {
    Sphere2D border;
    float density;
    std::vector<std::pair<int, float>> types;
    int patch_type_id;
    int biome_id;
    Patch(){};
    Patch(glm::vec2 pos, Biome::PatchDesc &patchDesc, int biome_id, int patch_type_id);
  };

  void generate_plants_cells(SceneGenerationContext &ctx, std::vector<int> cell_ids);
  void remove_trees_from_scene(SceneGenerationContext &ctx, std::vector<int> &ids);
  void prepare_tree_prototypes(SceneGenerationContext &ctx);
};

struct SceneGenerationContext
{
  SceneGenerationContext() : objects_bvh(false){};
  SceneGenerationContext &operator=(SceneGenerationContext &&ctx) = default;
  Scene scene;
  Block settings;
  BVH objects_bvh;
  BiomeMap biome_map;
  GroveMask global_mask;
  GrovePacker packer;
  std::vector<Cell> cells;
  std::vector<scene_gen::Patch> trees_patches;
  std::vector<scene_gen::Patch> grass_patches;
  int cells_x, cells_y;
  float hmap_pixel_size, biome_map_pixel_size;
  glm::vec2 heightmap_size, grass_field_size, full_size, cell_size, center, start_pos;
  glm::vec3 center3;

  int base_biome_id = 0;
  bool inited = false;
};
struct Cell
{
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
};