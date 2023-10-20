#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "render/world_renderer.h"
#include "generation/metainfo_manager.h"
#include "generation/scene_generator_helper.h"
#include <chrono>

Texture save_as_texture_RGBA8(const Field_2d &f, float mnv = 0, float mxv = 0)
{
      if (mnv == mxv)
        {
            mnv = f.min_val;
            mxv = f.max_val;
        }
        if (mnv <= mxv)
        {
            mnv = 0;
            mxv = 1;
        }

        if (!f.data)
            return engine::textureManager->empty();
        unsigned char *image_data = new unsigned char[4*(2*f.w + 1)*(2*f.h+1)];
        for (int i=-f.h;i<=f.h;i++)
        {
            for (int j=-f.w;j<=f.w;j++)
            {
                float val = (f.get(i,j) - mnv)/(mxv - mnv);
                int pos = (j+f.w)*(2*f.h + 1) + i + f.h;
                image_data[4*pos]   = 255*val;
                image_data[4*pos+1] = 255*val;
                image_data[4*pos+2] = 255*val;
                image_data[4*pos+3] = 255;
            }
        }

        Texture t = engine::textureManager->create_texture(2*f.w+1, 2*f.h+1, GL_RGBA8, 1, image_data, GL_RGBA);

        delete[] image_data;
        return t;
}

void RenderCmdExecutor::execute(int max_cmd_count)
{
    int cmd_left = max_cmd_count;
    while (!renderCmdBuffer->empty() && cmd_left != 0)
    {
        auto &cmd = renderCmdBuffer->front();
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        switch (cmd.type)
        {
        case RC_UPDATE_HMAP:
            if (genCtx.scene.heightmap)
                worldRenderer.set_heightmap(*(genCtx.scene.heightmap));
            else 
                worldRenderer.remove_heightmap();
            break;
        case RC_UPDATE_OBJECTS:
            worldRenderer.remove_all_instanced_models();
            worldRenderer.add_instanced_models(genCtx.scene.instanced_models);
            break;
        case RC_INIT_RENDER:
            if (!worldRenderer.is_inited())
            {
                worldRenderer.init(appCtx.windows_height, appCtx.window_width, cmd.args);
            }
            break;
        case RC_GLOBAL_PARAMS_UPDATE:
        {
            worldRenderer.debugInfo.set_bool("render_grid_debug", false);
            worldRenderer.debugInfo.set_bool("render_grove_mask_debug", false);
            worldRenderer.debugInfo.set_bool("render_biome_mask_debug", false);
            worldRenderer.debugInfo.set_bool("render_bvh_debug", false);
            worldRenderer.debugInfo.set_vec4("grid_params", glm::vec4(genCtx.start_pos, genCtx.cell_size));
            glm::vec2 gm_st_ps = glm::vec2(genCtx.global_mask.get_borders().x, genCtx.global_mask.get_borders().y);
            glm::vec2 gm_st_sz = glm::vec2((2*genCtx.global_mask.get_grid_size().x + 1) * genCtx.global_mask.get_cell_size(),
                                           (2*genCtx.global_mask.get_grid_size().y + 1) * genCtx.global_mask.get_cell_size());
            worldRenderer.debugInfo.set_vec4("debug_tex_scale", glm::vec4(gm_st_ps, 1.0f/gm_st_sz));
        }
            break;
        case RC_UPDATE_DEBUG_PARAMS:
            worldRenderer.debugInfo.add_detalization(cmd.args);
            break;
        case RC_UPDATE_TREES:
          if (worldRenderer.is_inited())
          {
            worldRenderer.set_grove(genCtx.scene.grove, AABB2D(genCtx.center - genCtx.grass_field_size, genCtx.center + genCtx.grass_field_size),
                                    metainfoManager->see_all_tree_types());
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
              Model *m = new Model();
              visualizer::visualize_light_voxels(genCtx.cells[cell_id].voxels_small, m);
              m->update();
              worldRenderer.debug_models.back().m  = m; 
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
                  {
                    Model *m = new Model();
                    visualizer::visualize_light_voxels(genCtx.cells[cell_id].voxels_small, m);
                    m->update();
                    dm.m = m;
                  }
                }
              }
            }
          }
          break;
        case RC_SAVE_GROVE_MASK_TO_TEXTURE:
          if (worldRenderer.debugInfo.get_bool("render_grove_mask_debug", false))
          {
            Texture t = save_as_texture_RGBA8(genCtx.global_mask);
            auto it = worldRenderer.debug_textures.find("grove_mask");
            if (it != worldRenderer.debug_textures.end())
            {
              engine::textureManager->delete_tex(it->second);
              it->second = t;
            }
            else
            {
              worldRenderer.debug_textures.emplace("grove_mask", t);
            }
          }
          break;
        case RC_SAVE_BIOME_MASK_TO_TEXTURE:
          if (worldRenderer.debugInfo.get_bool("render_biome_mask_debug", false))
          {
            Texture t = genCtx.biome_map.save_as_texture_RGBA8();
            auto it = worldRenderer.debug_textures.find("biome_mask");
            if (it != worldRenderer.debug_textures.end())
            {
              engine::textureManager->delete_tex(it->second);
              it->second = t;
            }
            else
            {
              worldRenderer.debug_textures.emplace("biome_mask", t);
            }
          }
          break;
        case RC_VISUALIZE_BVH_DEBUG:
          if (worldRenderer.debugInfo.get_bool("render_bvh_debug", false))
          {
            uint64_t model_id = SceneGenHelper::pack_id(0, (int)Scene::DEBUG_MODEL, 2, 0);
            for (auto &dm : worldRenderer.debug_models)
            {
              if (dm.id == model_id)
                delete dm.m;
            }
            std::vector<AABB> boxes;
            std::vector<glm::vec3> colors;
            int cnt = 0;
            auto func = [&](const std::pair<AABB, uint64_t> &p)
            {
              cnt = (cnt + 1)%27+1;
              boxes.push_back(p.first);
              colors.push_back(glm::vec3(0.5*(cnt % 3), 0.5*(cnt / 3 % 3), 0.5*(cnt / 9 % 3)));
            };
            genCtx.objects_bvh.iterate_over_intersected_bboxes(AABB(glm::vec3(-1e9, -1e9, -1e9),
                                                                    glm::vec3(1e9, 1e9, 1e9)),
                                                              func, false);
            worldRenderer.debug_models.emplace_back();
            worldRenderer.debug_models.back().apply_light = false;
            Model *m = new Model();
            visualizer::visualize_aabb(boxes, m, colors);
            m->update();
            worldRenderer.debug_models.back().m = m;
            worldRenderer.debug_models.back().id = model_id;
          }
          break;
        case RC_REMOVE_BVH_DEBUG:
        {
          uint64_t model_id = SceneGenHelper::pack_id(0, (int)Scene::DEBUG_MODEL, 2, 0);
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
        case RC_UPDATE_GRASS:
          worldRenderer.set_grass(genCtx.scene.grass);
          break;
        case RC_FINISH:
          appCtx.event.quit = true;
          break;
        default:
          logerr("RenderCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
          break;
        }
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        float ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        logerr("%s took %.3f ms", ToString(cmd.type), ms);
        renderCmdBuffer->pop();
        cmd_left--;
    }
    if (worldRenderer.is_inited())
        render();
}

void RenderCmdExecutor::render()
{
    worldRenderer.set_resolution(engine::view->WIDTH, engine::view->HEIGHT);
    worldRenderer.set_forced_LOD(appCtx.forced_LOD);
    worldRenderer.set_render_mode(appCtx.render_mode);
    RenderReadbackInputData rrid;
    RenderReadbackData rrd;
    rrid.cursor_screen_pos =  glm::vec2(appCtx.mousePos.x/engine::view->WIDTH, 1-(appCtx.mousePos.y)/engine::view->HEIGHT);
    worldRenderer.set_readback(rrid);
    worldRenderer.render(1, appCtx.camera);
    rrd = worldRenderer.get_readback();
    appCtx.mouseWorldPosType = rrd.cursor_world_pos_type;
}