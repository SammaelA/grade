#pragma once
#include "cmd_buffers.h"
#include "render/world_renderer.h"
#include "generation/scene_generator.h"
#include "app.h"

class InputCmdExecutor
{
public:
    void execute(int max_cmd_count = -1);
};

class GenerationCmdExecutor
{
public:
    void execute(int max_cmd_count = -1);
    GenerationCmdExecutor(AppContext &aCtx, SceneGenerator::SceneGenerationContext &gCtx):
    appCtx(aCtx),
    genCtx(gCtx)
    {

    }
private:
    AppContext &appCtx;
    SceneGenerator::SceneGenerationContext &genCtx;
};

class RenderCmdExecutor
{
public:
    RenderCmdExecutor(AppContext &aCtx, SceneGenerator::SceneGenerationContext &gCtx):
    appCtx(aCtx),
    genCtx(gCtx)
    {

    }
    void execute(int max_cmd_count = -1);
private:
    void render();
    WorldRenderer worldRenderer;
    AppContext &appCtx;
    SceneGenerator::SceneGenerationContext &genCtx;
};