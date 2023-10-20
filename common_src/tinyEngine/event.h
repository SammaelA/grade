#pragma once
#include "SDL2/SDL.h"
#include "view.h"

struct Scroll
{
  bool posx, posy, negx, negy;
  void reset()
  {
    posx = posy = negx = negy = false;
  }
};

class Event
{
public:
  Event(){};
  void input(); // Take inputs and add them to stack

  bool quit = false;
  bool fullscreenToggle = false;

  // Keyboard Events
  std::unordered_map<SDL_Keycode, bool> active;
  std::deque<SDL_Keycode> press;

  // Movement Events
  bool mousemove = false;
  SDL_MouseMotionEvent mouse;
  SDL_MouseWheelEvent mouseWheel;
  // Clicking Events
  std::unordered_map<Uint8, bool> click; // Active Buttons
  std::deque<Uint8> clicked;             // Button Events

  Scroll scroll;
  SDL_Event windowEvent; // Window Resizing Event
  bool windowEventTrigger = false;
};