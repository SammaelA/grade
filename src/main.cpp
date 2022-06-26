#include "cities_generator/cities_gen_main.h"
#include "tinyEngine/TinyEngine.h"
#include "cmd_buffers.h"
#include "app.h"
#include "graphics_utils/modeling.h"
#include "generation/scene_generator.h"
#include "cmd_executors.h"
#include "gui.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
CommandBuffer<InputCommands> inputCmdBuffer;
CommandBuffer<GenerationCommands> genCmdBuffer;
CommandBuffer<RenderCommands> renderCmdBuffer;

extern void read_from_console_nonblock();
extern void read_commands_from_string(std::string &block_str);

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
    Tiny::window("Procedural Tree", appContext.WIDTH, appContext.HEIGHT);
    BlkManager man;
    Block textures_list;
    man.load_block_from_file("resources.blk", textures_list);
    textureManager = TextureManager("./resources/textures/", textures_list);
    model_loader::load_default_blk();

    SceneGenerator::SceneGenerationContext sceneGenerationContext;
    sceneGenerationContext.scene = new Scene();

    InputCmdExecutor ice = InputCmdExecutor(sceneGenerationContext);
    GenerationCmdExecutor gce = GenerationCmdExecutor(sceneGenerationContext);
    RenderCmdExecutor rce = RenderCmdExecutor(appContext, sceneGenerationContext);
    GUI gui = GUI(appContext, sceneGenerationContext);
    InputHandler inputHandler = InputHandler(appContext, sceneGenerationContext);
    if (argc == 2)
    {
      std::string str(argv[1]);
      read_commands_from_string(str);
    }
    Tiny::event.handler = [&]()
    { 
      inputHandler.handle_input(Tiny::event);
    };
    Tiny::view.pipeline = [&]()
    {
        appContext.fpsCounter.tick();
        read_from_console_nonblock();
        ice.execute();
        gce.execute();
        rce.execute();
    };
    Tiny::view.interface = [&]() {
        gui.render_debug_settings();
        gui.render_cell_info();
    };

    Tiny::loop([&]() {});
    Tiny::quit();
}