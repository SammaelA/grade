#pragma once

enum InputCommands
{
    IC_GEN_HMAP,
    IC_ADD_OBJECT,
    IC_CLEAR_SCENE,
    IC_INIT_SCENE,
    IC_REMOVE_OBJECT,
    IC_COMMANDS_COUNT
};

enum GenerationCommands
{
    GC_GEN_HMAP,
    GC_ADD_OBJECT,
    GC_CLEAR_SCENE,
    GC_INIT_SCENE,
    GC_REMOVE_BY_ID,
    GC_COMMANDS_COUNT
};

enum RenderCommands
{
    RC_UPDATE_HMAP,
    RC_UPDATE_OBJECTS,
    RC_INIT_RENDER,
    RC_COMMANDS_COUNT
};