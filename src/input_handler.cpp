#include "input_handler.h"
#include "core/grove.h"
#include "generation/scene_generation.h"
#include "tinyEngine/camera.h"
#include "core/tree.h"
#include "app.h"
#include "cmd_executors.h"
#include "generation/scene_generator_helper.h"
#include "common_utils/matrix_transform.h"

void InputHandler::handle_input(Event &event)
{
  if(event.windowEventTrigger)
  {
    engine::view->WIDTH = event.windowEvent.window.data1;
    engine::view->HEIGHT = event.windowEvent.window.data2;
  }

  if (ctx.frames_from_last_input < 25) //do not use hotkeys when we are typing something
  {
    for (auto &e : event.active)
      event.active[e.first] = false;
    return;
  }
  
  ctx.window_width = engine::view->WIDTH;
  ctx.windows_height = engine::view->HEIGHT;
  
  float sensitivity = 0.1;
  float speed = 1;
  float mouse_scroll_speed = 5;

  float nx = event.mouse.x;
  float ny = event.mouse.y;
  float xoffset = nx - ctx.mousePos.x;
  float yoffset = ctx.mousePos.y - ny;
  ctx.mousePos = float2(nx, ny);
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
    float3 front;
    front.x = cos(LiteMath::to_radians(ctx.camera.yaw)) * cos(LiteMath::to_radians(ctx.camera.pitch));
    front.y = sin(LiteMath::to_radians(ctx.camera.pitch));
    front.z = sin(LiteMath::to_radians(ctx.camera.yaw)) * cos(LiteMath::to_radians(ctx.camera.pitch));
    ctx.camera.front = normalize(front);
  }

  if (event.active[SDLK_LALT])
    ctx.camera.pos += mouse_scroll_speed*(event.mouseWheel.y + 0.0f)*ctx.camera.front;

  //Pause Toggle
  float3 cameraPerp = normalize(cross(ctx.camera.front, ctx.camera.up));
  if (event.active[SDLK_t])
  {
    logerr("camera pos %f,%f,%f",ctx.camera.pos.x,ctx.camera.pos.y,ctx.camera.pos.z);
    logerr("camera front %f,%f,%f",ctx.camera.front.x,ctx.camera.front.y,ctx.camera.front.z);
    event.active[SDLK_t] = false;
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
    inputCmdBuffer->push(InputCommands::IC_GEN_HMAP);
    logerr("create hmap 1");
    event.active[SDLK_h] = false;
  }
  if (event.active[SDLK_ESCAPE])
  {
    inputCmdBuffer->push(InputCommands::IC_EXIT);
    event.active[SDLK_ESCAPE] = false;
  }
  if ((event.active[SDLK_LSHIFT] || event.active[SDLK_RSHIFT]) && event.active[SDLK_p])
  {
    Block b;
    b.add_vec4("world_pos_type", ctx.mouseWorldPosType);
    b.add_string("type_name", ctx.active_tree_type);
    bool can_plant_tree_here = (SceneGenHelper::is_terrain(ctx.mouseWorldPosType) && 
          abs(ctx.mouseWorldPosType.x - genCtx.start_pos.x) < genCtx.heightmap_size.x &&
          abs(ctx.mouseWorldPosType.z - genCtx.start_pos.y) < genCtx.heightmap_size.y);
    if (can_plant_tree_here)
      inputCmdBuffer->push(InputCommands::IC_PLANT_TREE_IMMEDIATE, b);
    event.active[SDLK_p] = false;
  }
  if (event.active[SDLK_o])
  {
    Block b;
    float4x4 transform = LiteMath::scale(LiteMath::translate(float4x4(), to_float3(ctx.mouseWorldPosType)), float3(ctx.cur_obj_scale));
    float4x4 transform2 = LiteMath::eulerAngleXYZ(ctx.cur_obj_angles.x, ctx.cur_obj_angles.y, ctx.cur_obj_angles.z);
    transform = transform*transform2;
    b.set_string("name", ctx.active_object_name);
    b.set_bool("on_terrain", ctx.cur_object_on_terrain);
    b.set_mat4("transform", transform);
    inputCmdBuffer->push(InputCommands::IC_ADD_OBJECT, b);
    event.active[SDLK_o] = false;
  }
  if (event.active[SDLK_b])
  {
    if (ctx.biome_brush >= 0)
      ctx.biome_brush = -1;
    event.active[SDLK_b] = false;
  }
  if (ctx.biome_brush >= 0)
  {
    if (event.click[SDL_BUTTON_LEFT] && SceneGenHelper::is_terrain(ctx.mouseWorldPosType))
    {
      Block b;
      b.set_vec3("pos", to_float3(ctx.mouseWorldPosType));
      b.set_double("outer_radius", ctx.biome_brush_size);
      b.set_double("inner_radius", 0.6*ctx.biome_brush_size);
      b.set_int("id", ctx.biome_brush);
      inputCmdBuffer->push(IC_SET_BIOME_ROUND, b);
      event.click[SDL_BUTTON_LEFT] = false;
    }
    else if (event.click[SDL_BUTTON_RIGHT] && SceneGenHelper::is_terrain(ctx.mouseWorldPosType))
    {
      Block b;
      b.set_vec3("pos", to_float3(ctx.mouseWorldPosType));
      b.set_double("outer_radius", ctx.biome_brush_size);
      b.set_double("inner_radius", 0.6*ctx.biome_brush_size);
      b.set_int("id", -1);
      inputCmdBuffer->push(IC_SET_BIOME_ROUND, b);
      event.click[SDL_BUTTON_RIGHT] = false;
    }
  }
  else
  {
    if (event.click[SDL_BUTTON_RIGHT])
    {
      float2 pos_xz = float2(ctx.mouseWorldPosType.x, ctx.mouseWorldPosType.z);
      int2 c_ij = to_int2((pos_xz - genCtx.start_pos)/genCtx.cell_size);
      int cell_id = c_ij.x*genCtx.cells_y + c_ij.y;
      if (cell_id >= 0 && cell_id < genCtx.cells.size())
      {
        if (cell_id == ctx.active_cell_id)
          ctx.active_cell_id = -1;
        else
        {
          auto it = ctx.cells.find(cell_id);
          if (it == ctx.cells.end())
            ctx.cells.emplace(cell_id, AppContext::CellUiInfo());
          ctx.active_cell_id = cell_id;
        }
      }
      event.click[SDL_BUTTON_RIGHT] = false;
    }
  }

  if(!event.press.empty())
  {
    if(event.press.back() == SDLK_F11)//Toggle fullscreen
    {   
      engine::view->fullscreen = !engine::view->fullscreen;
      if(!engine::view->fullscreen) 
        SDL_SetWindowFullscreen(engine::view->gWindow, 0);
      else 
        SDL_SetWindowFullscreen(engine::view->gWindow, SDL_WINDOW_FULLSCREEN_DESKTOP);
    }

    event.press.pop_back();
  }

  if(!event.clicked.empty()) event.clicked.pop_back();  //Reset Event Triggers
  event.scroll.reset();
  event.mousemove = false;
  event.windowEventTrigger = false;
  event.mouseWheel = SDL_MouseWheelEvent();
}