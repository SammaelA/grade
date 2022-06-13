#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "render/world_renderer.h"
#include "tinyEngine/TinyEngine.h"

void RenderCmdExecutor::execute(int max_cmd_count)
{
    Block default_render_settings;
    if (!worldRenderer.is_inited())
        worldRenderer.init(appCtx.HEIGHT, appCtx.WIDTH, default_render_settings);

    int cmd_left = max_cmd_count;
    while (!renderCmdBuffer.empty() && cmd_left != 0)
    {
        auto cmd = renderCmdBuffer.pop();
        switch (cmd.type)
        {
        case RC_UPDATE_HMAP:
            if (genCtx.scene->heightmap)
                worldRenderer.set_heightmap(*(genCtx.scene->heightmap));
            break;
        
        default:
            logerr("RenderCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
            break;
        }
        cmd_left--;
    }

    render();
}

void RenderCmdExecutor::render()
{
    worldRenderer.set_resolution(Tiny::view.WIDTH, Tiny::view.HEIGHT);
    worldRenderer.set_forced_LOD(appCtx.forced_LOD);
    worldRenderer.set_render_mode(appCtx.render_mode);
    worldRenderer.render(1, appCtx.camera);
}