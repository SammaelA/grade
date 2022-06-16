#include "gui.h"
#include "third_party/imgui/imgui.h"
#include "cmd_buffers.h"
#include "common_utils/utility.h"

void gui::render_debug_settings()
{
    constexpr int options_cnt = 1;
    static bool checks[options_cnt] = {false};
    static bool prev_checks[options_cnt] = {false};
    static std::string check_names[options_cnt] = {"Show grid"};
    static std::string options_names[options_cnt] = {"render_grid_debug"};
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