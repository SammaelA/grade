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
            renderCmdBuffer.push(RC_UPDATE_OBJECTS);
            break;
        case IC_CLEAR_SCENE:
            if (genCtx.inited)
            {
                genCmdBuffer.push(GC_CLEAR_SCENE);
                renderCmdBuffer.push(RC_UPDATE_HMAP);
                renderCmdBuffer.push(RC_UPDATE_OBJECTS);
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
            break;
        case IC_REMOVE_OBJECT:
        for (int i=0;i<cmd.args.size();i++)
            {
            uint64_t u_id = cmd.args.get_uint64(i, 0);
            BlkManager man;
            std::string s;
            man.save_block_to_string(s, cmd.args);
            logerr("%lu %s", u_id, s.c_str());
            if (u_id > 0)
            {
                unsigned a,b,c,d;
                SceneGenHelper::unpack_id(u_id, a,b,c,d);
                if (b == Scene::ObjCategories::SIMPLE_OBJECT)
                {
                    Block rm_ids;
                    //-1 in last field means that all objects of that type should be removed
                    rm_ids.add_ivec4("remove_mask", glm::ivec4(a,b,c,d));
                    logerr("remove by mask %d %d %d %d", a,b,c,d);
                    genCmdBuffer.push(GC_REMOVE_BY_ID, rm_ids);
                    renderCmdBuffer.push(RC_UPDATE_OBJECTS);
                }
            }
            }
            break;
        default:
            logerr("InputCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
            break;
        }
        inputCmdBuffer.pop();
        cmd_left--;
    }
}