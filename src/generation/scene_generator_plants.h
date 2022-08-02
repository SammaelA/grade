#pragma once
#include "core/scene.h"
#include "scene_generator.h"
namespace scene_gen
{
  void generate_plants_cells(SceneGenerationContext &ctx, std::vector<int> cell_ids);
  void remove_trees_from_scene(SceneGenerationContext &ctx, std::vector<int> &ids);
};