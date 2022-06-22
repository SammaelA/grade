#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "render/world_renderer.h"
#include "tinyEngine/TinyEngine.h"
#include "generation/metainfo_manager.h"
#include "generation/scene_generator_helper.h"

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
        {
            worldRenderer.debugInfo.set_bool("render_grid_debug", false);
            worldRenderer.debugInfo.set_bool("render_grove_mask_debug", false);
            worldRenderer.debugInfo.set_vec4("grid_params", glm::vec4(genCtx.start_pos, genCtx.cell_size));
            glm::vec2 gm_st_ps = glm::vec2(genCtx.global_mask->get_borders().x, genCtx.global_mask->get_borders().y);
            glm::vec2 gm_st_sz = glm::vec2((2*genCtx.global_mask->get_grid_size().x + 1) * genCtx.global_mask->get_cell_size(),
                                           (2*genCtx.global_mask->get_grid_size().y + 1) * genCtx.global_mask->get_cell_size());
            worldRenderer.debugInfo.set_vec4("debug_tex_scale", glm::vec4(gm_st_ps, 1.0f/gm_st_sz));
        }
            break;
        case RC_UPDATE_DEBUG_PARAMS:
            worldRenderer.debugInfo.add_detalization(cmd.args);
            break;
        case RC_UPDATE_TREES:
          if (worldRenderer.is_inited())
          {
            GroveGenerationData ggd;
            ggd.types = metainfoManager.get_all_tree_types();
            worldRenderer.set_grove(genCtx.scene->grove, ggd);
          }
          break;
        case RC_VISUALIZE_VOXELS_DEBUG:
        {
          int cell_id = cmd.args.get_int("cell_id",-1);
          if (cell_id >= 0)
          {
            logerr("visualize voxels in cell %d", cell_id);
            uint64_t model_id = SceneGenHelper::pack_id(0, (int)Scene::DEBUG_MODEL, 1, cell_id);
            bool found = false;
            for (auto &dm : worldRenderer.debug_models)
            {
              if (dm.id == model_id)
              {
                found = true;
                break;
              }
            }
            if (!found && genCtx.cells[cell_id].voxels_small)
            {
              worldRenderer.debug_models.emplace_back();
              worldRenderer.debug_models.back().apply_light = false;
              worldRenderer.debug_models.back().m = visualizer::visualize_light_voxels(genCtx.cells[cell_id].voxels_small);
              worldRenderer.debug_models.back().id = model_id;
            }
          }
          break;
        }
        case RC_REMOVE_VOXELS_DEBUG:
        {
          int cell_id = cmd.args.get_int("cell_id",-1);
          if (cell_id >= 0)
          {
            logerr("remove voxels in cell %d", cell_id);
            uint64_t model_id = SceneGenHelper::pack_id(0, (int)Scene::DEBUG_MODEL, 1, cell_id);
            for (auto dm = worldRenderer.debug_models.begin(); dm != worldRenderer.debug_models.end(); dm++)
            {
              if (dm->id == model_id)
              {
                delete dm->m;
                worldRenderer.debug_models.erase(dm);
                break;
              }
            }
          }
          break;
        }
        case RC_UPDATE_CELL:
          for (int i=0;i<cmd.args.size();i++)
          {
            int cell_id = cmd.args.get_int(i,-1);
            if (cell_id >= 0)
            {
              uint64_t model_id = SceneGenHelper::pack_id(0, (int)Scene::DEBUG_MODEL, 1, cell_id);
              for (auto &dm : worldRenderer.debug_models)
              {
                if (dm.id == model_id)
                {
                  delete dm.m;
                  if (genCtx.cells[cell_id].voxels_small)
                    dm.m = visualizer::visualize_light_voxels(genCtx.cells[cell_id].voxels_small);
                }
              }
            }
          }
          break;
        case RC_SAVE_GROVE_MASK_TO_TEXTURE:
          if (worldRenderer.debugInfo.get_bool("render_grove_mask_debug", false))
          {
            Texture t = genCtx.global_mask->save_as_texture_RGBA8();
            auto it = worldRenderer.debug_textures.find("grove_mask");
            if (it != worldRenderer.debug_textures.end())
            {
              textureManager.delete_tex(it->second);
              it->second = t;
            }
            else
            {
              worldRenderer.debug_textures.emplace("grove_mask", t);
            }
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
    RenderReadbackInputData rrid;
    RenderReadbackData rrd;
    rrid.cursor_screen_pos =  glm::vec2(appCtx.mousePos.x/Tiny::view.WIDTH, 1-(appCtx.mousePos.y)/Tiny::view.HEIGHT);
    worldRenderer.set_readback(rrid);
    worldRenderer.render(1, appCtx.camera);
    rrd = worldRenderer.get_readback();
    appCtx.mouseWorldPosType = rrd.cursor_world_pos_type;
}