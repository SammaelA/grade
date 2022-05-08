#include "app.h"

FpsCounter::FpsCounter()
{
    t1 = std::chrono::steady_clock::now();
    t_prev = std::chrono::steady_clock::now();
}
void FpsCounter::tick()
{
    t_prev = t1;
    t1 = std::chrono::steady_clock::now();
    float frame_time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t_prev).count();
    frame_time = MAX(frame_time,0.1);
    average_fps = mu*average_fps + (1 - mu)*(1000/frame_time);
    frame++;
}

std::function<void(AppContext &, Event &)> eventHandler = [](AppContext &ctx, Event &event)
{
  float nx = event.mouse.x;
  float ny = event.mouse.y;
  GLfloat xoffset = nx - ctx.mousePos.x;
  GLfloat yoffset = ctx.mousePos.y - ny;

  GLfloat sensitivity = 0.1;
  xoffset *= sensitivity;
  yoffset *= sensitivity;

  ctx.camera.yaw += xoffset;
  ctx.camera.pitch += yoffset;

  if (ctx.camera.pitch > 89.0f)
    ctx.camera.pitch = 89.0f;
  if (ctx.camera.pitch < -89.0f)
    ctx.camera.pitch = -89.0f;
  glm::vec3 front;
  front.x = cos(glm::radians(ctx.camera.yaw)) * cos(glm::radians(ctx.camera.pitch));
  front.y = sin(glm::radians(ctx.camera.pitch));
  front.z = sin(glm::radians(ctx.camera.yaw)) * cos(glm::radians(ctx.camera.pitch));
  ctx.camera.front = glm::normalize(front);
  ctx.mousePos = glm::vec2(nx, ny);
  //Pause Toggle
  float speed = 0.75;
  glm::vec3 cameraPerp = glm::normalize(glm::cross(ctx.camera.front, ctx.camera.up));
  if (event.active[SDLK_1])
  {
    ctx.camera = Camera();
    ctx.camera.pos = glm::vec3(-84,61.35,6);
    ctx.camera.front = glm::vec3(0.88,-0.05,0.47);
  }
  if (event.active[SDLK_2])
  {
    ctx.camera = Camera();
    ctx.camera.pos = glm::vec3(-223.5,68,23);
    ctx.camera.front = glm::vec3(0,70,0) - ctx.camera.pos;
  }
  if (event.active[SDLK_3])
  {
    ctx.camera = Camera();
    ctx.camera.pos = glm::vec3(-168,67.8,2.1);
    ctx.camera.front = glm::vec3(0.996,-0.056,-0.07);
  }
  if (event.active[SDLK_4])
  {
    ctx.camera = Camera();
    ctx.camera.pos = glm::vec3(0,70,-400);
    ctx.camera.front = glm::vec3(0,70,0) - ctx.camera.pos;
  }
  if (event.active[SDLK_5])
  {
    logerr("camera pos %f,%f,%f",ctx.camera.pos.x,ctx.camera.pos.y,ctx.camera.pos.z);
    logerr("camera front %f,%f,%f",ctx.camera.front.x,ctx.camera.front.y,ctx.camera.front.z);
    event.active[SDLK_5] = false;
  }
  if (event.active[SDLK_w])
    ctx.camera.pos += speed * ctx.camera.front;
  if (event.active[SDLK_s])
    ctx.camera.pos -= speed * ctx.camera.front;

  if (event.active[SDLK_a])
    ctx.camera.pos -= speed * cameraPerp;
  if (event.active[SDLK_d])
    ctx.camera.pos += speed * cameraPerp;

  if (event.active[SDLK_e])
    ctx.camera.pos += speed * ctx.camera.up;
  if (event.active[SDLK_c])
    ctx.camera.pos -= speed * ctx.camera.up;
  if (event.active[SDLK_y])
    ctx.forced_LOD = -1;
  if (event.active[SDLK_u])
    ctx.forced_LOD = 0;
  if (event.active[SDLK_i])
    ctx.forced_LOD = 1;
  if (event.active[SDLK_o])
    ctx.forced_LOD = 2;
  if (event.active[SDLK_p])
    ctx.forced_LOD = 3;
  if (event.active[SDLK_LEFTBRACKET])
    ctx.forced_LOD = 4;
  if (event.active[SDLK_RIGHTBRACKET])
    ctx.forced_LOD = 5;

  if (event.active[SDLK_m])
  {
    if (ctx.benchmark_grove_current == -1)
    {
      ctx.render_mode++;
      if (ctx.render_mode > ctx.MAX_RENDER_MODE)
        ctx.render_mode = ctx.ARRAY_TEX_DEBUG_RENDER_MODE;
      logerr("render mode %d",ctx.render_mode);
    }
    else
    {
      ctx.benchmark_grove_needed = ctx.benchmark_grove_current + 1;
    }
    event.active[SDLK_m] = false;
  }
  if (event.active[SDLK_n])
  {
    if (ctx.render_mode <= ctx.DEBUG_RENDER_MODE)
      ctx.debug_tex++;
    event.active[SDLK_n] = false;
  }
  if (event.active[SDLK_b])
  {
    if (ctx.groveRendererDebugParams.need_focus_model)
      ctx.groveRendererDebugParams.model_focused++;
    if (ctx.render_mode == ctx.ARRAY_TEX_DEBUG_RENDER_MODE)
      ctx.debug_layer++;
    event.active[SDLK_b] = false;
  }
    if (event.active[SDLK_v])
  {
    if (ctx.groveRendererDebugParams.need_focus_model)
      ctx.groveRendererDebugParams.model_focused--;
    if (ctx.render_mode == ctx.ARRAY_TEX_DEBUG_RENDER_MODE)
      ctx.debug_layer--;
    event.active[SDLK_v] = false;
  }
  if (event.active[SDLK_f])
  {
    ctx.groveRendererDebugParams.need_focus_model = !ctx.groveRendererDebugParams.need_focus_model;
    ctx.render_mode = -1;
    event.active[SDLK_f] = false;
  }
  if (event.active[SDLK_r])
  {
    ctx.regeneration_needed = true;
    event.active[SDLK_r] = false;
  }
  if (event.active[SDLK_t])
  {
    ctx.add_generation_needed = true;
    event.active[SDLK_t] = false;
  }
  if (event.active[SDLK_ESCAPE])
  {
    exit(0);
    event.active[SDLK_ESCAPE] = false;
  }
  if (!event.press.empty())
  {

  }
};