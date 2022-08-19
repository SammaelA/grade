#pragma once
#include <functional>
#include <initializer_list>

#include "../third_party/imgui/imgui.h"                    //Interface Dependencies
#include "../third_party/imgui/imgui_impl_sdl.h"
#include "../third_party/imgui/imgui_impl_opengl3.h"

#include "tinyEngine/resources.h"                                //Rendering Dependencies
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_mixer.h>
#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"

#include <sstream>                                  //File / Console IO
#include <iostream>
#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "texture.h"              //Utility Classes
#include "shader.h"
#include "model.h"

#include "view.h"
#include "event.h"
#include "audio.h"
#include "common_utils/utility.h"
#include <chrono>

/* TINY ENGINE */

namespace Tiny {

extern View view;           //Window and Interface  (Requires Initialization)
extern Event event;         //Event Handler
extern Audio audio;         //Audio Processor       (Requires Initialization)

bool window(std::string windowName, int width, int height);
void quit();

template<typename F, typename... Args>
void loop(F function, Args&&... args){
  while(!event.quit){    
    if(Tiny::audio.enabled) audio.process();      //Audio Processor

    function(args...);      //User-defined Game Loop

    if(Tiny::view.enabled)  view.render();        //Render View
    if(Tiny::view.enabled){
      event.input();        //Get Input
      event.handle(view);   //Call the event-handling system
    }
    else
    {
      logerr("view not enabled\n");
    }
  }
}

}
