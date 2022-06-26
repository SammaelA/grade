#include "app.h"
#include "cmd_buffers.h"

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

void InputHandler::handle_input(Event &event)
{
  float sensitivity = 0.1;
  float speed = 1;
  float mouse_scroll_speed = 5;

  float nx = event.mouse.x;
  float ny = event.mouse.y;
  float xoffset = nx - ctx.mousePos.x;
  float yoffset = ctx.mousePos.y - ny;
  ctx.mousePos = glm::vec2(nx, ny);
  if (event.active[SDLK_LALT] || event.click[SDL_BUTTON_MIDDLE] || ctx.free_camera)
  {
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
  }

  if (event.active[SDLK_LALT])
    ctx.camera.pos += mouse_scroll_speed*(event.mouseWheel.y + 0.0f)*ctx.camera.front;

  //Pause Toggle
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
  //if (event.active[SDLK_p])
  //  ctx.forced_LOD = 3;
  if (event.active[SDLK_LEFTBRACKET])
    ctx.forced_LOD = 4;
  if (event.active[SDLK_RIGHTBRACKET])
    ctx.forced_LOD = 5;
  if (event.active[SDLK_h])
  {
    //ctx.save_to_hydra = true;
    inputCmdBuffer.push(InputCommands::IC_GEN_HMAP);
    logerr("create hmap 1");
    event.active[SDLK_h] = false;
  }
  if (event.active[SDLK_ESCAPE])
  {
    exit(0);
    event.active[SDLK_ESCAPE] = false;
  }
  if (event.active[SDLK_p])
  {
    Block b;
    b.add_vec4("world_pos_type", ctx.mouseWorldPosType);
    b.add_string("type_name", ctx.active_tree_type);
    inputCmdBuffer.push(InputCommands::IC_PLANT_TREE_IMMEDIATE, b);
    event.active[SDLK_p] = false;
  }

  if (event.click[SDL_BUTTON_RIGHT])
  {
    glm::vec2 pos_xz = glm::vec2(ctx.mouseWorldPosType.x, ctx.mouseWorldPosType.z);
    glm::ivec2 c_ij = (pos_xz - genCtx.start_pos)/genCtx.cell_size;
    int cell_id = c_ij.x*genCtx.cells_y + c_ij.y;
    if (cell_id >= 0 && cell_id < genCtx.cells.size())
    {
      if (cell_id == ctx.active_cell_id)
        ctx.active_cell_id = -1;
      else
      {
        auto it = ctx.cells.find(cell_id);
        if (it == ctx.cells.end())
          ctx.cells.emplace(cell_id, AppContext::CellUiInfo()).first;
        ctx.active_cell_id = cell_id;
      }
    }
    event.click[SDL_BUTTON_RIGHT] = false;
  }
  if (!event.press.empty())
  {

  }
};