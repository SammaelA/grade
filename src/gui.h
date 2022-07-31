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
  void render_model_creation_info();
  void text_input();

  void read_commands_from_string(std::string &block_str);
  void read_from_console_nonblock();
private:
  AppContext &appCtx;
  const SceneGenerator::SceneGenerationContext &genCtx;
};