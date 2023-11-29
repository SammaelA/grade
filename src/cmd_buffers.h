#pragma once
#include <typeinfo>
#include "commands.h"
#include "common_utils/blk.h"
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
        commands.front().args.set_int("cmd_code", (int)commands.front().type);
        commands.front().args.set_string("cmd_name", ToString(commands.front().type));
        cmd_log.add_block("cmd_"+std::to_string(commands.front().id), &(commands.front().args));
        commands.pop();
        current_command_id++;
        save_block_to_file(std::string("../logs/") + typeid(T).name() + std::string("_log.blk"), cmd_log);
    }
    CommandBuffer() = default;
    ~CommandBuffer()= default;
private:
    std::queue<Command> commands;
    Block cmd_log;
    unsigned current_command_id = 0;
    unsigned next_command_id = 0;
};

extern CommandBuffer<InputCommands> *inputCmdBuffer;
extern CommandBuffer<GenerationCommands> *genCmdBuffer;
extern CommandBuffer<RenderCommands> *renderCmdBuffer;