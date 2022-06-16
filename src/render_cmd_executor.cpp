#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "render/world_renderer.h"
#include "tinyEngine/TinyEngine.h"

void RenderCmdExecutor::execute(int max_cmd_count)
{
    int cmd_left = max_cmd_count;
    while (!renderCmdBuffer.empty() && cmd_left != 0)
    {
        auto &cmd = renderCmdBuffer.front();
        switch (cmd.type)
        {
        case RC_UPDATE_HMAP:
            if (genCtx.scene->heightmap)
                worldRenderer.set_heightmap(*(genCtx.scene->heightmap));
            else 
                worldRenderer.remove_heightmap();
            break;
        case RC_UPDATE_OBJECTS:
            worldRenderer.remove_all_instanced_models();
            worldRenderer.add_instanced_models(genCtx.scene->instanced_models);
            break;
        case RC_INIT_RENDER:
            if (!worldRenderer.is_inited())
            {
                Block default_render_settings;
                worldRenderer.init(appCtx.HEIGHT, appCtx.WIDTH, default_render_settings);
            }
            break;
        case RC_GLOBAL_PARAMS_UPDATE:
            worldRenderer.debugInfo.set_bool("render_grid_debug", false);
            worldRenderer.debugInfo.set_vec4("grid_params", glm::vec4(genCtx.center, genCtx.cell_size));
            break;
        case RC_UPDATE_DEBUG_PARAMS:
            worldRenderer.debugInfo.add_detalization(cmd.args);
            break;
        default:
            logerr("RenderCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
            break;
        }
        renderCmdBuffer.pop();
        cmd_left--;
    }
    if (worldRenderer.is_inited())
        render();
}

void RenderCmdExecutor::render()
{
    worldRenderer.set_resolution(Tiny::view.WIDTH, Tiny::view.HEIGHT);
    worldRenderer.set_forced_LOD(appCtx.forced_LOD);
    worldRenderer.set_render_mode(appCtx.render_mode);
    RenderReadbackInputData rrid;
    RenderReadbackData rrd;
    rrid.cursor_screen_pos =  glm::vec2(appCtx.mousePos.x/Tiny::view.WIDTH, 1-(appCtx.mousePos.y)/Tiny::view.HEIGHT);
    worldRenderer.set_readback(rrid);
    worldRenderer.render(1, appCtx.camera);
    rrd = worldRenderer.get_readback();
    float dist = rrd.cursor_on_geometry ? length(appCtx.camera.pos - rrd.cursor_world_pos) : -1;
    appCtx.mouseWorldPosDist = glm::vec4(rrd.cursor_world_pos, dist);
}