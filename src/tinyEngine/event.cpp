#include "event.h"
#include "common_utils/utility.h"
void Event::input()
{
  if (SDL_PollEvent(&in) == 0)
    return;
  ImGui_ImplSDL2_ProcessEvent(&in);

  switch (in.type)
  {
  case SDL_QUIT:
    quit = true;
    break;
  case SDL_KEYDOWN:
    active[in.key.keysym.sym] = true;
    break;
  case SDL_KEYUP:
    active[in.key.keysym.sym] = false;
    press.push_front(in.key.keysym.sym);
    break;
  case SDL_MOUSEWHEEL:
    scroll.posy = (in.wheel.y > 0.99);
    scroll.negy = (in.wheel.y < -0.99);
    scroll.posx = (in.wheel.x > 0.99);
    scroll.negx = (in.wheel.x < -0.99);
    mouseWheel = in.wheel;
    break;
  case SDL_MOUSEMOTION:
    mouse = in.motion;
    mousemove = true;
    break;
  case SDL_MOUSEBUTTONDOWN:
    click[in.button.button] = true;
    break;
  case SDL_MOUSEBUTTONUP:
    click[in.button.button] = false;
    clicked.push_front(in.button.button);
    break;
  case SDL_WINDOWEVENT:
    if (in.window.event == SDL_WINDOWEVENT_RESIZED)
    {
      windowEvent = in;
      windowEventTrigger = true;
    }
    break;
  default:
    break;
  }
}
