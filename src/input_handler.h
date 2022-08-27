#pragma once
#include "core/grove.h"
#include "generation/scene_generation.h"
#include "tinyEngine/camera.h"
#include "core/tree.h"
#include "app.h"

class InputHandler
{
public:
  InputHandler(AppContext &app_ctx, const SceneGenerationContext &gen_ctx):
  ctx(app_ctx),
  genCtx(gen_ctx)
  {

  };
  void handle_input(Event &e);
private:
  AppContext &ctx;
  const SceneGenerationContext &genCtx;
};