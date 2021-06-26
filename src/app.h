#pragma once
#include <functional>
#include <chrono>
#include "grove.h"
#include "tinyEngine/camera.h"
#include "tree.h"
#include "tinyEngine/event.h"
class FpsCounter
{
    float average_fps;
    uint64_t frame = 0;
    float mu = 0.99;
    std::chrono::steady_clock::time_point t1, t_prev;
public:
    FpsCounter();
    void tick();
    float get_average_fps() {return average_fps;}
    int get_frame_n() {return frame;}
};
struct AppContext
{
    const int WIDTH = 1200;
    const int HEIGHT = 800;
    const int DEBUG_RENDER_MODE = -2;
    const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
    const int MAX_RENDER_MODE = 2;
    const float fov = glm::radians(90.0f);
    int forced_LOD = -1;
    int render_mode = -1;
    int debug_tex = 0;
    int debug_layer = 0;
    glm::vec2 mousePos = glm::vec2(-1, -1);
    glm::mat4 projection = glm::perspective(fov, (float)WIDTH / HEIGHT, 1.0f, 3000.0f);

    GroveRendererDebugParams groveRendererDebugParams;
    FpsCounter fpsCounter;
    Camera camera;
    DirectedLight light;
};

extern std::function<void(AppContext &, Event &)> eventHandler;