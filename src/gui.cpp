#include "gui.h"
#include "third_party/imgui/imgui.h"
#include "cmd_buffers.h"
#include "common_utils/utility.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/modeling.h"
#include "third_party/icons.h"
#include <map>

GUI::GUI(AppContext &app_ctx, const SceneGenerator::SceneGenerationContext &gen_ctx) : appCtx(app_ctx),
                                                                                       genCtx(gen_ctx)
{
  ImGuiIO &io = ImGui::GetIO();
  io.Fonts->AddFontDefault();

  // merge in icons from Font Awesome
  static const ImWchar icons_ranges[] = {ICON_MIN_FA, ICON_MAX_16_FA, 0};
  ImFontConfig icons_config;
  icons_config.MergeMode = true;
  icons_config.PixelSnapH = true;
  io.Fonts->AddFontFromFileTTF(FONT_ICON_FILE_NAME_FAS, 16.0f, &icons_config, icons_ranges);
};

void GUI::render_main_toolbar()
{
  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoMove;
  window_flags |= ImGuiWindowFlags_NoResize;
  window_flags |= ImGuiWindowFlags_NoCollapse;
  window_flags |= ImGuiWindowFlags_NoTitleBar;
  window_flags |= ImGuiWindowFlags_NoScrollbar;
  ImGui::SetNextWindowPos(ImVec2(0, -10), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImVec2(appCtx.window_width, 60), ImGuiCond_Always);
  ImGui::SetNextWindowContentSize(ImVec2(appCtx.window_width, 60));
  ImGui::Begin(" ", nullptr, window_flags);
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(15, 15));
  ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(10, 10));

  bool press_1 = ImGui::Button(ICON_FA_GEAR);ImGui::SameLine();
  bool press_2 = ImGui::Button(ICON_FA_HOUSE);ImGui::SameLine();
  bool press_3 = ImGui::Button(ICON_FA_TREE);ImGui::SameLine();
  bool press_4 = ImGui::Button(ICON_FA_COMPUTER);

  ImGui::PopStyleVar(2);
  ImGui::End();
  {
    static bool show = false;
    if (press_1)
      show = !show;
    if (show)
      render_debug_settings();
  }
  {
    static bool show = false;
    if (press_2)
      show = !show;
    if (show)
      render_model_creation_info();
  }
  {
    static bool show = false;
    if (press_3)
      show = !show;
    if (show)
      render_tree_plant_info();
  }
  {
    static bool show = false;
    if (press_4)
      show = !show;
    if (show)
      text_input();
  }

  if (appCtx.active_cell_id >= 0)
    render_cell_info();
}

void GUI::render_debug_settings()
{
    constexpr int options_cnt = 3;
    static bool checks[options_cnt] = {false, false, false};
    static bool prev_checks[options_cnt] = {false, false, false};
    static std::string check_names[options_cnt] = {"Show grid", "Show global mask", "Show bvh"};
    static std::string options_names[options_cnt] = {"render_grid_debug", "render_grove_mask_debug", "render_bvh_debug"};
    ImGui::Begin("Render debug params");
    for (int i=0;i<options_cnt;i++)
        ImGui::Checkbox(check_names[i].c_str(), &(checks[i]));
    ImGui::End();

    bool change = false;
    Block b;
    for (int i=0;i<options_cnt;i++)
    {
        if (checks[i] != prev_checks[i])
        {
            change = true;
            b.add_bool(options_names[i], checks[i]);
        }
        prev_checks[i] = checks[i];
    }

    if (change)
    {
        inputCmdBuffer.push(InputCommands::IC_UPDATE_RENDER_DEBUG_PARAMS, b);
    }
}

void GUI::render_cell_info()
{
  if (appCtx.active_cell_id < 0)
    return;
  auto &cell = appCtx.cells.at(appCtx.active_cell_id);
  bool visualize_prev = cell.visualize_voxels_array;
  ImGui::Begin("Cell info"); 
  ImGui::Text("Cell id: %d", appCtx.active_cell_id);
  ImGui::Checkbox("Show voxels array", &(cell.visualize_voxels_array));
  bool remove_trees = ImGui::Button("Remove trees");
  ImGui::End();

  if (visualize_prev && !cell.visualize_voxels_array)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer.push(InputCommands::IC_REMOVE_VOXELS_DEBUG, b);
  }
  else if (!visualize_prev && cell.visualize_voxels_array)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer.push(InputCommands::IC_VISUALIZE_VOXELS_DEBUG, b);    
  }
  if (remove_trees)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer.push(InputCommands::IC_CLEAR_CELL, b);  
  }
}

void GUI::render_model_creation_info()
{
  static int cur_item = -1;
  int prev_cur_item = cur_item;
  static std::vector<std::string> names;
  static const char *types[256];
  static int types_cnt = 0;
  static bool update = true;

  static float angle_x = 0, angle_y = 0, angle_z = 0;
  static float scale = 1;
  static bool on_terrain = false;

  ImGui::Begin("Create object"); 
  ImGui::ListBox("Available models", &cur_item, types, types_cnt);
  ImGui::SliderAngle("Rotate X", &angle_x);
  ImGui::SliderAngle("Rotate Y", &angle_y);
  ImGui::SliderAngle("Rotate Z", &angle_z);
  ImGui::SliderFloat("Scale", &scale, 1, 1000);
  ImGui::Checkbox("On terrain", &on_terrain);
  appCtx.cur_obj_angles = glm::vec3(angle_x, angle_y, angle_z);
  appCtx.cur_obj_scale = scale;
  appCtx.cur_object_on_terrain = on_terrain;
  bool b_update = ImGui::Button("Update objects list");
  ImGui::End();

  if ((update || b_update) && model_loader::obj_models_blk_loaded)
  {
    names.clear();
    for (int i=0;i<256;i++)
      types[i] = nullptr;
    types_cnt = model_loader::obj_models_blk.size();
    names.reserve(types_cnt);
    for (int j=0;j<types_cnt;j++)
    {
      if (j>255)
        break;
      names.push_back(model_loader::obj_models_blk.get_name(j));
      types[j] = names.back().c_str();
    }
    update = false;
  }

  if (cur_item != prev_cur_item)
  {
    prev_cur_item = cur_item;
    if (cur_item >= 0)
      appCtx.active_object_name = types[cur_item];
  }
}

void GUI::render_tree_plant_info()
{
  static int cur_item = -1;
  int prev_cur_item = cur_item;
  const char * types[256] = {nullptr};
  auto &types_map = metainfoManager.see_all_tree_type_names();
  int i=0;
  for (auto &p : types_map)
  {
    if (i>255)
      break;
    types[i] = p.first.c_str();
    i++;
  }
  ImGui::Begin("Plant tree"); 
  ImGui::ListBox("Available Types", &cur_item, types, types_map.size());
  ImGui::End();
  if (cur_item != prev_cur_item)
  {
    prev_cur_item = cur_item;
    if (cur_item >= 0)
      appCtx.active_tree_type = types[cur_item];
  }
}