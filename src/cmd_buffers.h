#pragma once
#include "commands.h"
#include "save_utils/blk.h"
#include <queue>

template <typename T>
class CommandBuffer
{
public:
    struct Command
    {
        T type;
        unsigned long id;
        Block args;
        Command(T _type, unsigned long _id)
        {
            type = _type;
            id = _id;
        }
        Command(T _type, unsigned long _id, Block &_args)
        {
            type = _type;
            id = _id;
            args = Block();
            args.copy(&_args);
        }
    };
    void push(T cmd_type)
    {
        commands.emplace(cmd_type, next_command_id);
        next_command_id++;
    }
    void push(T cmd_type, Block &args)
    {
        commands.emplace(cmd_type, next_command_id, args);
        next_command_id++;
    }
    bool empty()
    {
        return commands.empty();
    }
    Command &front()
    {
        return commands.front();
    }
    void pop()
    {
        current_command_id++;
        commands.pop();
    }
private:
    std::queue<Command> commands;
    int current_command_id = -1;
    int next_command_id = 0;
};

extern CommandBuffer<InputCommands> inputCmdBuffer;
extern CommandBuffer<GenerationCommands> genCmdBuffer;
extern CommandBuffer<RenderCommands> renderCmdBuffer;