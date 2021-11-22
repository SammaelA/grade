#pragma once
#include "grove_generator.h"
#include "core/scene.h"

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
    LightVoxelsCube *create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                                         AABB &influence_box);
    AABB get_influence_AABB(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                            Heightmap &h);
    void generate_grove();
    SceneGenerationContext &ctx;
};

struct Cell
{
  enum CellStatus
  {
    EMPTY,
    WAITING,
    BORDER,
    FINISHED
  };
  GrovePrototype prototype;
  LightVoxelsCube *voxels_small = nullptr;
  int id = -1;
  AABB2D bbox;
  CellStatus status;
  std::vector<int> depends;//list of waiting cell (ids) that will use voxels from this cell
  std::vector<int> depends_from;
  std::vector<int> grass_patches;//set and used by grass generator
  AABB influence_bbox;
  explicit Cell(CellStatus _status = CellStatus::EMPTY) {status = _status;}
};