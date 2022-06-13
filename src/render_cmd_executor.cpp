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
    worldRenderer.render(1, appCtx.camera);
}