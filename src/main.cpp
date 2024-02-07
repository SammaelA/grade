#include "cmd_buffers.h"
#include "app.h"
#include "tree_utils/tree_modeling.h"
#include "generation/scene_generation.h"
#include "cmd_executors.h"
#include "gui.h"
#include "generation/metainfo_manager.h"
#include "input_handler.h"
#include <csignal>

MetainfoManager *metainfoManager = nullptr;
CommandBuffer<InputCommands> *inputCmdBuffer = nullptr;
CommandBuffer<GenerationCommands> *genCmdBuffer = nullptr;
CommandBuffer<RenderCommands> *renderCmdBuffer = nullptr;

void defaultSignalHandler(int signum)
{
  logerr("Interrupt signal received. Closing the program");
  exit(signum);  
}

int main(int argc, char *argv[])
{
    signal(SIGINT, defaultSignalHandler);  
    signal(SIGTSTP, defaultSignalHandler); 

    AppContext appContext;
    View view;
    view.lineWidth = 1.0f;
    view.init("Procedural Tree", appContext.window_width, appContext.windows_height);
    view.hide_window();
    engine::view = &view;

    model_loader::load_default_blk();
    Block textures_list;
    load_block_from_file("resources.blk", textures_list);
    TextureManager textureManager = TextureManager("./resources/textures/", textures_list);
    CommandBuffer<InputCommands> input_cmd_buffer;
    CommandBuffer<GenerationCommands> gen_cmd_buffer;
    CommandBuffer<RenderCommands> render_cmd_buffer;
    MetainfoManager metainfo_manager;
    engine::textureManager = &textureManager;
    inputCmdBuffer = &input_cmd_buffer;
    genCmdBuffer = &gen_cmd_buffer;
    renderCmdBuffer = &render_cmd_buffer;
    metainfoManager = &metainfo_manager;

    SceneGenerationContext sceneGenerationContext;
    InputCmdExecutor       ice(sceneGenerationContext);
    GenerationCmdExecutor  gce(sceneGenerationContext);
    RenderCmdExecutor      rce(appContext, sceneGenerationContext);
    GUI                    gui(appContext, sceneGenerationContext);
    InputHandler           inputHandler(appContext, sceneGenerationContext);

    if (argc >= 2)
    {
      std::string str(argv[1]);
      gui.read_commands_from_string(str);
    }
    while(!appContext.event.quit)
    {     
      view.audio.process();
      appContext.event.input();
      inputHandler.handle_input(appContext.event);
      
      appContext.fpsCounter.tick();
      ice.execute();
      gce.execute();
      rce.execute();
      gui.render();

     view.next_frame();
    }
    view.quit();
    logerr("exiting normally");
}