#include "cities_generator/cities_gen_main.h"
#include "tinyEngine/TinyEngine.h"
#include "cmd_buffers.h"
#include "app.h"
#include "graphics_utils/modeling.h"
#include "generation/scene_generation.h"
#include "cmd_executors.h"
#include "gui.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
CommandBuffer<InputCommands> inputCmdBuffer;
CommandBuffer<GenerationCommands> genCmdBuffer;
CommandBuffer<RenderCommands> renderCmdBuffer;

int main(int argc, char *argv[])
{
    //if (argc == 1)
    //    return cities_generator_main();
    //else
    //    return parser::parser_main(argc, argv);
    //return old_main(argc,argv);
    AppContext appContext;

    glewInit();
    Tiny::view.lineWidth = 1.0f;
    Tiny::window("Procedural Tree", appContext.window_width, appContext.windows_height);
    
    Block textures_list;
    load_block_from_file("resources.blk", textures_list);
    textureManager = TextureManager("./resources/textures/", textures_list);
    model_loader::load_default_blk();

    SceneGenerationContext sceneGenerationContext;
    InputCmdExecutor ice = InputCmdExecutor(sceneGenerationContext);
    GenerationCmdExecutor gce = GenerationCmdExecutor(sceneGenerationContext);
    RenderCmdExecutor rce = RenderCmdExecutor(appContext, sceneGenerationContext);
    GUI gui = GUI(appContext, sceneGenerationContext);
    InputHandler inputHandler = InputHandler(appContext, sceneGenerationContext);
    if (argc >= 2)
    {
      std::string str(argv[1]);
      gui.read_commands_from_string(str);
    }

    Tiny::view.interface = [&]() {
        gui.render_main_toolbar();
        gui.read_from_console_nonblock();
    };

    while(!Tiny::event.quit)
    {    
      if(Tiny::audio.enabled) 
        Tiny::audio.process();

      if(Tiny::view.enabled)
      {
        appContext.fpsCounter.tick();
        ice.execute();
        gce.execute();
        rce.execute();
        Tiny::view.drawInterface();
        SDL_GL_SwapWindow(Tiny::view.gWindow);
        Tiny::event.input();
        inputHandler.handle_input(Tiny::event);
      }
    }
    Tiny::quit();
}