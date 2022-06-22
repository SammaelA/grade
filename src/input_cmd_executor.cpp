#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "generation/scene_generator_helper.h"

void InputCmdExecutor::execute(int max_cmd_count)
{
  int cmd_left = max_cmd_count;
  while (!inputCmdBuffer.empty() && cmd_left != 0)
  {
    auto &cmd = inputCmdBuffer.front();
    switch (cmd.type)
    {
    case IC_GEN_HMAP:
      genCmdBuffer.push(GC_GEN_HMAP, cmd.args);
      renderCmdBuffer.push(RC_UPDATE_HMAP);
      break;
    case IC_ADD_OBJECT:
      genCmdBuffer.push(GC_ADD_OBJECT, cmd.args);
      genCmdBuffer.push(GC_UPDATE_GLOBAL_MASK);
      renderCmdBuffer.push(RC_UPDATE_OBJECTS);
      renderCmdBuffer.push(RC_VISUALIZE_BVH_DEBUG);
      renderCmdBuffer.push(RC_SAVE_GROVE_MASK_TO_TEXTURE);
      break;
    case IC_CLEAR_SCENE:
      if (genCtx.inited)
      {
        genCmdBuffer.push(GC_CLEAR_SCENE);
        genCmdBuffer.push(GC_UPDATE_GLOBAL_MASK);
        renderCmdBuffer.push(RC_UPDATE_HMAP);
        renderCmdBuffer.push(RC_UPDATE_OBJECTS);
        renderCmdBuffer.push(RC_SAVE_GROVE_MASK_TO_TEXTURE);
      }
      break;
    case IC_INIT_SCENE:
      renderCmdBuffer.push(RC_INIT_RENDER);
      if (genCtx.inited)
      {
        genCmdBuffer.push(GC_CLEAR_SCENE);
        renderCmdBuffer.push(RC_UPDATE_HMAP);
        renderCmdBuffer.push(RC_UPDATE_OBJECTS);
      }
      genCmdBuffer.push(GC_INIT_SCENE, cmd.args);
      renderCmdBuffer.push(RC_GLOBAL_PARAMS_UPDATE);
      break;
    case IC_REMOVE_OBJECT:
      for (int i = 0; i < cmd.args.size(); i++)
      {
        uint64_t u_id = cmd.args.get_uint64(i, 0);
        BlkManager man;
        std::string s;
        man.save_block_to_string(s, cmd.args);
        if (u_id > 0)
        {
          unsigned a, b, c, d;
          SceneGenHelper::unpack_id(u_id, a, b, c, d);
          if (b == Scene::ObjCategories::SIMPLE_OBJECT)
          {
            Block rm_ids;
            //-1 in last field means that all objects of that type should be removed
            rm_ids.add_ivec4("remove_mask", glm::ivec4(a, b, c, d));
            genCmdBuffer.push(GC_REMOVE_BY_ID, rm_ids);
            genCmdBuffer.push(GC_UPDATE_GLOBAL_MASK);
            renderCmdBuffer.push(RC_UPDATE_OBJECTS);
            renderCmdBuffer.push(RC_VISUALIZE_BVH_DEBUG);
            renderCmdBuffer.push(RC_SAVE_GROVE_MASK_TO_TEXTURE);
          }
        }
      }
      break;
    case IC_UPDATE_RENDER_DEBUG_PARAMS:
      renderCmdBuffer.push(RC_UPDATE_DEBUG_PARAMS, cmd.args);
      if (cmd.args.get_bool("render_grove_mask_debug", false))
      {
        renderCmdBuffer.push(RC_SAVE_GROVE_MASK_TO_TEXTURE);
      }
      if (cmd.args.get_id("render_bvh_debug") >= 0)
      {
        if (cmd.args.get_bool("render_bvh_debug", false))
          renderCmdBuffer.push(RC_VISUALIZE_BVH_DEBUG);
        else 
          renderCmdBuffer.push(RC_REMOVE_BVH_DEBUG);
      }
      break;
    case IC_PLANT_TREE:
    case IC_PLANT_TREE_IMMEDIATE:
    {
      if (!genCtx.inited)
        break;
      glm::vec4 world_pos_type = cmd.args.get_vec4("world_pos_type");
      if (SceneGenHelper::is_terrain(world_pos_type) && 
          abs(world_pos_type.x - genCtx.start_pos.x) < genCtx.heightmap_size.x &&
          abs(world_pos_type.z - genCtx.start_pos.y) < genCtx.heightmap_size.y)
      {
        //we can plant the tree actually
        std::string tree_name = "small_oak_simplified";//TODO - choose tree name
        Block b;
        glm::vec2 pos = glm::vec2(world_pos_type.x, world_pos_type.z);
        b.add_vec2("pos", pos);
        b.add_string("type_name", tree_name);
        genCmdBuffer.push(GC_PLANT_TREE, b);

        if (cmd.type == IC_PLANT_TREE_IMMEDIATE)
        {
          glm::ivec2 c_ij = (pos - genCtx.start_pos) / genCtx.cell_size;
          Block cb;
          cb.add_int("cell_id", c_ij.x*genCtx.cells_y + c_ij.y);
          genCmdBuffer.push(GC_GEN_TREES_CELL, cb);
          renderCmdBuffer.push(RC_UPDATE_CELL, cb);
          renderCmdBuffer.push(RC_UPDATE_TREES);
        }
      }
      break;
    }
    case IC_GEN_ALL_PLANTED_TREES:
      if (genCtx.inited)
      {
        int cells_to_update = 0;
        Block cb = Block();
        for (int i=0;i<genCtx.cells_x;i++)
        {
          for (int j=0;j<genCtx.cells_y;j++)
          {
            const Cell &c = genCtx.cells[i*genCtx.cells_y + j];
            for (const auto &p : c.prototypes)
            {
              if (p.trees_count > 0 || p.preplanted_trees.size() > 0)
              {
                cb.add_int("cell_id", i*genCtx.cells_y + j);
                cells_to_update++;
              }
            }
          }
        }
        if (cells_to_update > 0)
        {
          genCmdBuffer.push(GC_GEN_TREES_CELL, cb);
          renderCmdBuffer.push(RC_UPDATE_CELL, cb);
          renderCmdBuffer.push(RC_UPDATE_TREES);
        }
      }
      break;
    case IC_VISUALIZE_VOXELS_DEBUG:
      renderCmdBuffer.push(RC_VISUALIZE_VOXELS_DEBUG, cmd.args);
      break;
    case IC_REMOVE_VOXELS_DEBUG:
      renderCmdBuffer.push(RC_REMOVE_VOXELS_DEBUG, cmd.args);
      break;
    default:
      logerr("InputCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
      break;
    }
    inputCmdBuffer.pop();
    cmd_left--;
  }
}