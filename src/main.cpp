#include "parser.h"
#include "cities_generator/cities_gen_main.h"
#include "tinyEngine/TinyEngine.h"
#include "cmd_buffers.h"
#include "app.h"
#include "graphics_utils/modeling.h"
#include "generation/scene_generator.h"
#include "cmd_executors.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
CommandBuffer<InputCommands> inputCmdBuffer;
CommandBuffer<GenerationCommands> genCmdBuffer;
CommandBuffer<RenderCommands> renderCmdBuffer;

extern void read_from_console_nonblock();

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
    Tiny::event.handler = [&]()
    { eventHandler(appContext, Tiny::event); };
    BlkManager man;
    Block textures_list;
    man.load_block_from_file("resources.blk", textures_list);
    textureManager = TextureManager("./resources/textures/", textures_list);
    ModelLoader::load_default_blk();

    SceneGenerator::SceneGenerationContext sceneGenerationContext;
    sceneGenerationContext.scene = new Scene();

    InputCmdExecutor ice = InputCmdExecutor(sceneGenerationContext);
    GenerationCmdExecutor gce = GenerationCmdExecutor(sceneGenerationContext);
    RenderCmdExecutor rce = RenderCmdExecutor(appContext, sceneGenerationContext);

    Tiny::view.pipeline = [&]()
    {
        read_from_console_nonblock();
        ice.execute();
        gce.execute();
        rce.execute();
    };
    Tiny::view.interface = [&]() {
        ImGui::Begin("Demo window");
        ImGui::Button("Hello!");
        ImGui::End();
    };

    Tiny::loop([&]() {});
    {};

    Tiny::quit();
}