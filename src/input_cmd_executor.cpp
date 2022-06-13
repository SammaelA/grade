#include "cmd_executors.h"
#include "common_utils/utility.h"

void InputCmdExecutor::execute(int max_cmd_count)
{
    int cmd_left = max_cmd_count;
    while (!inputCmdBuffer.empty() && cmd_left != 0)
    {
        auto cmd = inputCmdBuffer.pop();
        switch (cmd.type)
        {
        case IC_GEN_HMAP:
            genCmdBuffer.push(GC_GEN_HMAP, cmd.args);
            renderCmdBuffer.push(RC_UPDATE_HMAP);
            break;
        
        default:
            logerr("InputCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
            break;
        }
        cmd_left--;
    }
}