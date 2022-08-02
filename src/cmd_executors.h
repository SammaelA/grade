#pragma once
#include "cmd_buffers.h"
#include "render/world_renderer.h"
#include "generation/scene_generator.h"
#include "app.h"

class InputCmdExecutor
{
public:
    void execute(int max_cmd_count = -1);
    InputCmdExecutor(const SceneGenerationContext &gCtx):
    genCtx(gCtx)
    {

    }
private:
    const SceneGenerationContext &genCtx;
};

class GenerationCmdExecutor
{
public:
    void execute(int max_cmd_count = -1);
    GenerationCmdExecutor(SceneGenerationContext &gCtx):
    genCtx(gCtx)
    {

    }
private:
    SceneGenerationContext &genCtx;
};

class RenderCmdExecutor
{
public:
    RenderCmdExecutor(AppContext &aCtx, const SceneGenerationContext &gCtx):
    appCtx(aCtx),
    genCtx(gCtx)
    {

    }
    void execute(int max_cmd_count = -1);
private:
    void render();
    WorldRenderer worldRenderer;
    AppContext &appCtx;
    const SceneGenerationContext &genCtx;
};