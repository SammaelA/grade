#pragma once
#include <functional>
#include <chrono>
#include "core/grove.h"
#include "generation/scene_generator.h"
#include "tinyEngine/camera.h"
#include "core/tree.h"
#include "common_utils/sun.h"
#include "tinyEngine/event.h"
class FpsCounter
{
  float average_fps;
  uint64_t frame = 0;
  float mu = 0.99;
  std::chrono::steady_clock::time_point t1, t_prev;

public:
  FpsCounter();
  void tick();
  float get_average_fps() { return average_fps; }
  int get_frame_n() { return frame; }
};

struct AppContext
{
  const int WIDTH = 1000;
  const int HEIGHT = 1000;
  const int DEBUG_RENDER_MODE = -2;
  const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
  const int MAX_RENDER_MODE = 2;
  const float fov = glm::radians(90.0f);
  int forced_LOD = 4;
  int render_mode = 0;
  glm::vec2 mousePos = glm::vec2(-1, -1);
  glm::vec4 mouseWorldPosType = glm::vec4(0, 0, 0, -1); //-1 means that mouse is not on scene geometry
  FpsCounter fpsCounter;
  Camera camera;
  bool free_camera = false;

  struct CellUiInfo
  {
    bool visualize_voxels_array = false;
  };
  std::map<int, CellUiInfo> cells;
  int active_cell_id = -1;

  std::string active_tree_type = "small_oak_simplified";
};

class InputHandler
{
public:
  InputHandler(AppContext &app_ctx, const SceneGenerator::SceneGenerationContext &gen_ctx):
  ctx(app_ctx),
  genCtx(gen_ctx)
  {

  };
  void handle_input(Event &e);
private:
  AppContext &ctx;
  const SceneGenerator::SceneGenerationContext &genCtx;
};