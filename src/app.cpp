#include "app.h"
#include "cmd_buffers.h"
#include "generation/scene_generator_helper.h"
#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/euler_angles.hpp>

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
    bool can_plant_tree_here = (SceneGenHelper::is_terrain(ctx.mouseWorldPosType) && 
          abs(ctx.mouseWorldPosType.x - genCtx.start_pos.x) < genCtx.heightmap_size.x &&
          abs(ctx.mouseWorldPosType.z - genCtx.start_pos.y) < genCtx.heightmap_size.y);
    if (can_plant_tree_here)
      inputCmdBuffer.push(InputCommands::IC_PLANT_TREE_IMMEDIATE, b);
    event.active[SDLK_p] = false;
  }
  if (event.active[SDLK_o])
  {
    Block b;
    glm::mat4 transform = glm::scale(glm::translate(glm::mat4(1.0f), glm::vec3(ctx.mouseWorldPosType)), glm::vec3(ctx.cur_obj_scale));
    glm::mat4 transform2 = glm::eulerAngleXYZ(ctx.cur_obj_angles.x, ctx.cur_obj_angles.y, ctx.cur_obj_angles.z);
    transform = transform*transform2;
    b.set_string("name", ctx.active_object_name);
    b.set_bool("on_terrain", ctx.cur_object_on_terrain);
    b.set_mat4("transform", transform);
    inputCmdBuffer.push(InputCommands::IC_ADD_OBJECT, b);
    event.active[SDLK_o] = false;
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