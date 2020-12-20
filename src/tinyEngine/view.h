#pragma once
#include <deque>
#include <unordered_map>
#include <string>
#include "SDL2/SDL.h"
#include <functional>
#include <glm/glm.hpp>
#include "imgui/imgui.h"                    //Interface Dependencies
#include "imgui/imgui_impl_sdl.h"
#include "imgui/imgui_impl_opengl3.h"
#include <GL/glew.h>    

class View{
  using Handle = std::function<void()>;
  public:
    bool init(std::string windowName, int width, int height);
    void quit();
    bool enabled = false;

    unsigned int WIDTH, HEIGHT;

    SDL_Window* gWindow;        //Window Pointer
    SDL_GLContext gContext;     //Render Context

    ImGuiIO io;
    Handle interface = [](){};  //User defined Interface
    bool showInterface = false;
    void drawInterface();

    Handle pipeline = [](){};           //User defined Pipeline
    void render();
    void target(glm::vec3 clearcolor);  //Target main window for drawing

    bool fullscreen = false;    //Settings
    bool vsync = true;
    bool ccw = true;
    unsigned int antialias = 16;
    float lineWidth = 1.0f;
};