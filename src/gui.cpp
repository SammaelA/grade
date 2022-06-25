#include "gui.h"
#include "third_party/imgui/imgui.h"
#include "cmd_buffers.h"
#include "common_utils/utility.h"
#include <map>

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
}