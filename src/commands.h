#pragma once
#include "common_utils/utility.h"

DEFINE_ENUM_WITH_STRING_CONVERSIONS(InputCommands,
  (IC_GEN_HMAP)
  (IC_ADD_OBJECT)
  (IC_CLEAR_SCENE)
  (IC_INIT_SCENE)
  (IC_REMOVE_OBJECT)
  (IC_UPDATE_RENDER_DEBUG_PARAMS)
  (IC_PLANT_TREE)
  (IC_PLANT_TREE_IMMEDIATE)
  (IC_GEN_ALL_PLANTED_TREES)
  (IC_VISUALIZE_VOXELS_DEBUG)
  (IC_REMOVE_VOXELS_DEBUG)
  (IC_CLEAR_CELL)
  (IC_EXIT)
  (IC_TREE_GEN_PARAMETER_SELECTION)
  (IC_GEN_BIOME_MAP)
  (IC_SET_DEFAULT_BIOME)
  (IC_SET_BIOME_ROUND)
  (IC_PREPARE_PLANT_PROTOTYPES)
  (IC_REMOVE_ALL_PLANTS)
  (IC_EXPORT_SCENE_TO_HYDRA)
  (IC_SAVE_SCENE)
  (IC_LOAD_SCENE)
  (IC_INIT_RENDER)
  (IC_VISUALIZE_TREE_HYDRA)
  (IC_COMMANDS_COUNT)
)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(GenerationCommands,
  (GC_GEN_HMAP)
  (GC_ADD_OBJECT)
  (GC_CLEAR_SCENE)
  (GC_INIT_SCENE)
  (GC_REMOVE_BY_ID)
  (GC_PLANT_TREE)
  (GC_GEN_TREES_CELL)
  (GC_UPDATE_GLOBAL_MASK)
  (GC_REMOVE_TREES)
  (GC_TREE_GEN_PARAMETER_SELECTION)
  (GC_GEN_BIOME_MAP)
  (GC_SET_DEFAULT_BIOME)
  (GC_SET_BIOME_ROUND)
  (GC_PREPARE_PLANT_PROTOTYPES)
  (GC_REMOVE_GRASS_IN_CELLS)
  (GC_SAVE_SCENE)
  (GC_LOAD_SCENE)
  (GC_VISUALIZE_TREE_HYDRA)
  (GC_COMMANDS_COUNT)
)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(RenderCommands,
  (RC_UPDATE_HMAP)
  (RC_UPDATE_OBJECTS)
  (RC_INIT_RENDER)
  (RC_GLOBAL_PARAMS_UPDATE)
  (RC_UPDATE_DEBUG_PARAMS)
  (RC_UPDATE_TREES)
  (RC_VISUALIZE_VOXELS_DEBUG)
  (RC_REMOVE_VOXELS_DEBUG)
  (RC_UPDATE_CELL)
  (RC_SAVE_GROVE_MASK_TO_TEXTURE)
  (RC_VISUALIZE_BVH_DEBUG)
  (RC_REMOVE_BVH_DEBUG)
  (RC_SAVE_BIOME_MASK_TO_TEXTURE)
  (RC_UPDATE_GRASS)
  (RC_FINISH)
  (RC_COMMANDS_COUNT)
)