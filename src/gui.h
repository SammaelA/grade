#pragma once
#include "app.h"
#include "generation/scene_generator.h"
class GUI
{
public:
  GUI(AppContext &app_ctx, const SceneGenerator::SceneGenerationContext &gen_ctx):
  appCtx(app_ctx),
  genCtx(gen_ctx)
  {

  };
  void render_debug_settings();
  void render_cell_info();
  void render_tree_plant_info();
private:
  AppContext &appCtx;
  const SceneGenerator::SceneGenerationContext &genCtx;
};