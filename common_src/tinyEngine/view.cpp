#include "view.h"
#include "configuration.h"
#include "third_party/imgui/imgui.h"                    //Interface Dependencies
#include "third_party/imgui/imgui_impl_sdl.h"
#include "third_party/imgui/imgui_impl_opengl3.h"

#include "resources.h"                                //Rendering Dependencies
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_mixer.h>
#include "common_utils/LiteMath_ext.h"
#include "common_utils/matrix_transform.h"

#include <sstream>                                  //File / Console IO
#include <iostream>
#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include "engine.h"

View *engine::view = nullptr;

void View::hide_window()
{
  if (!window_hidden)
    SDL_HideWindow(gWindow);
  window_hidden = true;
}

void View::show_window()
{
  if (window_hidden)
    SDL_ShowWindow(gWindow);
  window_hidden = false;
}

bool View::init(std::string _name, int W, int H)
{
  glewInit();
  
  if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
      //Если тут возникает ошибшка, попробуйте запуск в том же терминале
      // комманды $ export SDL_VIDEO_X11_VISUALID= 
      //Или сразу отключайте мультисэмплинг
      printf("SDL could not initialize! Error: %s\n", SDL_GetError());
      return false;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG))
    {
      printf("SDL_Image could not initialize! Error: %s\n", IMG_GetError());
      return false;
    }

    if (TTF_Init() == -1)
    { // for some reason this is -1
      printf("SDL_ttf could not initialize! Error: %s\n", TTF_GetError());
      return false;
    }

  enabled = true;
  WIDTH = W;
  HEIGHT = H;

  // Core OpenGL Profile!
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

  //В случае segfault на базовых openGL коммандах попробуйте отключить мультисамплинг
  if (get_antialiasing())
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, get_antialiasing());

  // Initialize the Window and Context
  gWindow = SDL_CreateWindow(_name.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL);
  if (gWindow == NULL)
  {
    printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
    return false;
  }
  SDL_SetWindowResizable(gWindow, SDL_TRUE);
  gContext = SDL_GL_CreateContext(gWindow);

  SDL_GL_SetSwapInterval(vsync);
  glewExperimental = GL_TRUE; // Launch GLEW
  glewInit();

  IMGUI_CHECKVERSION(); // Setup ImGUI
  ImGui::CreateContext();
  io = ImGui::GetIO();
  (void)io;
  ImGui_ImplSDL2_InitForOpenGL(gWindow, gContext);
  ImGui_ImplOpenGL3_Init("#version 330 core");
  ImGui::StyleColorsCustom();

  if (get_antialiasing())
    glEnable(GL_MULTISAMPLE);
  glEnable(GL_DEPTH_TEST); // Setup Global OpenGL State!
  glDepthFunc(GL_LEQUAL);
  // glEnable(GL_ALPHA_TEST);
  glDisable(GL_CULL_FACE);
  if (ccw)
    glFrontFace(GL_CCW);
  else
    glFrontFace(GL_CW);
  glLineWidth(lineWidth);
  glEnable(GL_PROGRAM_POINT_SIZE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  if (!audio.init())
    { // Start the Audio Interface
      std::cout << "Failed to launch audio interface." << std::endl;
      return false;
    }

  return true;
}

void View::quit()
{
  audio.quit();
  ImGui_ImplOpenGL3_Shutdown(); // Shutdown ImGUI
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(gContext);
  SDL_DestroyWindow(gWindow);
    
  TTF_Quit();
  IMG_Quit();
  SDL_Quit();
}

void View::target(glm::vec3 clearcolor)
{
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glViewport(0, 0, WIDTH, HEIGHT);
  glClearColor(clearcolor.x, clearcolor.y, clearcolor.z, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void View::next_frame()
{
  SDL_GL_SwapWindow(gWindow);
}

unsigned int View::get_antialiasing()
{
  return VIEW_ANTIALIASING;
}