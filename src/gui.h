#pragma once
#include "app.h"
#include "generation/scene_generator.h"
class GUI
{
public:
  GUI(AppContext &app_ctx, const SceneGenerationContext &gen_ctx);
  void render_main_toolbar();
  void render_debug_settings();
  void render_cell_info();
  void render_tree_plant_info();
  void render_model_creation_info();
  void render_biome_toolbar();
  void text_input();
  void render_parameter_selection_menu();
  void blk_modification_interface(Block *b, const std::string &title);
  void read_commands_from_string(std::string &block_str);
  void read_from_console_nonblock();
  bool has_active_input() {return input_active;}
private:
  bool input_text(const char *label, std::string &text);
  bool input_active = false;
  AppContext &appCtx;
  const SceneGenerationContext &genCtx;
};