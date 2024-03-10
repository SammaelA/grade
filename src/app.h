#pragma once
#include <functional>
#include <chrono>
#include "core/grove.h"
#include "generation/scene_generation.h"
#include "tinyEngine/camera.h"
#include "core/tree.h"
#include "tinyEngine/event.h"
#include "tinyEngine/audio.h"
#include "tinyEngine/view.h"
#include "tinyEngine/engine.h"

struct AppContext
{
  Event event; //Event Handler

  const int DEBUG_RENDER_MODE = -2;
  const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
  const int MAX_RENDER_MODE = 2;
  const float fov = LiteMath::to_radians(90.0f);
  int forced_LOD = 4;
  int render_mode = 0;
  int window_width = 1000;
  int windows_height = 1000;

  glm::vec2 mousePos = glm::vec2(-1, -1);
  glm::vec4 mouseWorldPosType = glm::vec4(0, 0, 0, -1); //-1 means that mouse is not on scene geometry
  FpsCounter fpsCounter;
  Camera camera;
  bool free_camera = false;
  int frames_from_last_input = 1000;

  struct CellUiInfo
  {
    bool visualize_voxels_array = false;
  };
  std::map<int, CellUiInfo> cells;
  int active_cell_id = -1;

  std::string active_tree_type = "small_oak_simplified";
  std::string active_object_name = "farm_1";
  glm::vec3 cur_obj_angles = glm::vec3(0,0,0);
  bool cur_object_on_terrain = false;
  float cur_obj_scale = 1;

  int biome_brush = -1;
  float biome_brush_size = 40;
};
