#pragma once
#include "grove_generator.h"
#include "core/scene.h"
#include <mutex>
#include <atomic>

struct Cell;
class SceneGenerator
{
public:
    struct SceneGenerationContext
    {
        Scene *scene;
        Block settings;
        std::map<std::string,TreeTypeData> tree_types;
        std::map<std::string, GrassType> grass_types;
        GroveGenerationData global_ggd;
    };
    
    SceneGenerator(SceneGenerationContext &_ctx): ctx(_ctx){};
    void create_scene_auto();
private:
    //LightVoxelsCube *create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
    //                                     AABB &influence_box);
    //LightVoxelsCube *create_cell_small_voxels(Cell &c);
    //AABB get_influence_AABB(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
    //                        Heightmap &h);
    void generate_grove();
    SceneGenerationContext &ctx;
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
  GrovePrototype prototype;
  LightVoxelsCube *voxels_small = nullptr;//used for plants gen is dependant cells
  Field_2d *planar_occlusion = nullptr;//used for grass gen in this cell
  int id = -1;
  AABB2D bbox;
  std::atomic<CellStatus> status;
  std::vector<int> depends;//list of waiting cell (ids) that will use voxels from this cell
  std::vector<int> depends_from;
  std::vector<int> grass_patches;//set and used by grass generator
  AABB influence_bbox;
  explicit Cell(CellStatus _status = CellStatus::EMPTY) {status = _status;}
  Cell(const Cell &cell)
  {
    id = cell.id;
    prototype = cell.prototype;
    bbox = cell.bbox;
    status.store(cell.status);
    influence_bbox = cell.influence_bbox;
    depends = cell.depends;
    depends_from = cell.depends_from;

    if (voxels_small || planar_occlusion || !grass_patches.empty())
    {
      logerr("only copying of empty cells is allowed");
    }
  };
};