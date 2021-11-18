#include "parser.h"
#define GLEW_EXPERIMENTAL
#include "tinyEngine/TinyEngine.h"
#include "glm/trigonometric.hpp"
#include "tinyEngine/image.h"
//#include "tinyEngine/color.h"
//#include "tinyEngine/helper.h"
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include <string>
#include "tinyEngine/camera.h"
#include "graphics_utils/modeling.h"
#include "graphics_utils/texture_manager.h"
#include "common_utils/utility.h"
#include "core/grove.h"
#include "render/grove_renderer.h"
#include "save_utils/config.h"
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include "graphics_utils/terrain.h"
#include "render/shadow.h"
#include "generation/grove_packer.h"
#include "app.h"
#include "render/grass_renderer.h"
#include "tinyEngine/deffered_target.h"
#include "render/ambient_occlusion.h"
#include "tinyEngine/cubemap.h"
#include "common_utils/python_interaction.h"
#include "gltf_utils/general_gltf_writer.h"
#include <thread>
#include "generation/parameter_selection.h"
#include "save_utils/blk.h"
#include "clustering/clustering_benchmark.h"
#include "tree_generators/load_tree_structure.h"
#include "clustering/clustering_debug_status.h"
#include "generation/grove_generator.h"
#include "tree_generators/GE_generator.h"
#include "tree_generators/python_tree_gen.h"
#include "tree_generators/weber_penn_parameters.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/proctree.h"
#include "tree_generators/generated_tree.h"
#include "render/billboard_cloud_renderer.h"
#include "render/visualizer.h"
#include "render/world_renderer.h"
#include "generation/scene_generator.h"

namespace parser
{
//View Tiny::view;   //Window and Interface  (Requires Initialization)
//Event Tiny::event; //Event Handler
//Audio Tiny::audio; //Audio Processor       (Requires Initialization)
AppContext appContext;

bool generation_needed = false;
bool saving_needed = false;
bool loading_needed = false;
bool print_perf = false;
bool only_gen = false;
bool visualize_voxels = false;
bool visualize_initial_voxels = false;
bool statistics_run = false;
bool gltf_export = false;
bool need_initialization = true;
bool parameter_selection = false;
bool clustering_benchmark = false;

std::string generator_name = "default";
std::string generator_fixed_preset_name = "";
std::string parameter_selector_name = "default_selection";
std::string clustering_benchmark_path = "benchmark.blk";
struct StatRunLaunchParams
{
  int trees = 1;
  float max_ind_dist = 0.7;
} statRunLaunchParams;
GroveRenderer::Precision pres = GroveRenderer::MEDIUM;
std::string grove_type_name = "default";
std::string save_path = ".";
std::string load_path = ".";

int parse_arguments(int argc, char *argv[])
{
  int k = 1;
  while (k < argc)
  {
    if (std::string(argv[k]) == "-g")
    {
      generation_needed = true;
      if (argc == k + 1)
      {
        logerr("write grove type name after -g");
        return 1;
      }
      else
        grove_type_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-perf")
    {
      print_perf = true;
      k++;
    }
    else if (std::string(argv[k]) == "-parameter_selection")
    {
      parameter_selection = true;
      if (argc > k+1)
      {
        parameter_selector_name = std::string(argv[k+1]);
        k++;
      }
      k++;
    }
    else if (std::string(argv[k]) == "-benchmark")
    {
      clustering_benchmark = true;
      if (argc > k+1)
      {
        clustering_benchmark_path = std::string(argv[k+1]);
        k++;
      }
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_voxels")
    {
      visualize_voxels = true;
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_clusters")
    {
      clusteringDebugInfo.visualize_clusters = true;
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_initial_voxels")
    {
      visualize_initial_voxels = true;
      k++;
    }
    else if (std::string(argv[k]) == "-only_gen")
    {
      only_gen = true;
      k++;
    }
    else if (std::string(argv[k]) == "-export_gltf")
    {
      gltf_export = true;
      k++;
    }
    else if (std::string(argv[k]) == "-s")
    {
      saving_needed = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -s");
        return 1;
      }
      else
        save_path = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-l")
    {
      loading_needed = true;
      if (argc == k + 1)
      {
        logerr("write path to load from after -l");
        return 1;
      }
      else
        load_path = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-low_precision")
    {
      pres = GroveRenderer::LOW;
      k++;
    }
    else if (std::string(argv[k]) == "-high_precision")
    {
      pres = GroveRenderer::DEBUG;
      k++;
    }
    else if (std::string(argv[k]) == "-h")
    {
      logerr("-g <grove type name> -s <save path> -l <load path>");
      return 1;
    }
    else if (std::string(argv[k]) == "-stat_run")
    {
      if (argc <= k + 2)
      {
        logerr("-stat_run <num_tree> <max_individual distance>");
      }
      else
      {
        statRunLaunchParams.trees = atoi(argv[k + 1]);
        statRunLaunchParams.max_ind_dist = atof(argv[k + 2]);
        statistics_run = true;
        //only_gen = true;
        generation_needed = true;
        debug("stat run %d trees %f MID\n", statRunLaunchParams.trees, statRunLaunchParams.max_ind_dist);
        k += 3;
      }
    }
    else if (std::string(argv[k]) == "-generator")
    {
      bool ok = argc >= k+2;
      if (ok)
      {
        if (std::string(argv[k+1]) == "=")
          generator_name = std::string(argv[k+2]);
        else
          ok = false;
      }
      if (!ok)
      {
        logerr("-generator = <generator_name>");
        logerr("possible names : default, proctree, simple, python_tree_gen");
      }
      k += 3;
      if (generator_name == "python_tree_gen")
      {
        bool ok = argc >= k+2;
        if (ok && std::string(argv[k]) == "-parameters")
        {
          if (std::string(argv[k+1]) == "=")
            generator_fixed_preset_name = std::string(argv[k+2]);
          else
            ok = false;
        }
        else 
          ok = false;

        if (ok)
          k += 3;
      }
    }
    else if (std::string(argv[k]) == "-prepare_dataset")
    {
      clusteringDebugInfo.prepare_dataset = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -prepare_dataset");
        return 1;
      }
      else
        clusteringDebugInfo.dataset_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-save_csv")
    {
      clusteringDebugInfo.save_csv = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -save_csv");
        return 1;
      }
      else
        clusteringDebugInfo.csv_file_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-progress_bar")
    {
      show_progress = 1;
      k++;
    }
    else
    {
      logerr("unknown command \"%s\". Write -h for help", argv[k]);
      return 1;
    }
  }
  if (saving_needed && !generation_needed)
  {
    logerr("Error: -s should be used only with -g flag");
    return 1;
  }
  if (!generation_needed && !loading_needed)
  {
    need_initialization = false;
    //logerr("You choosed to do nothing. Program will be closed");
    //return 0;
  }
  return -1;
}
void load_tree_types(std::map<std::string,TreeTypeData> &tree_types)
{
  int id = 0;
  BlkManager man;
  Block ge_gen_types, my_gen_types;
  man.load_block_from_file("ge_gen_presets.blk",ge_gen_types);
  man.load_block_from_file("my_gen_presets.blk",my_gen_types);

  for (int i=0;i<ge_gen_types.size();i++)
  {
    Block *bl = ge_gen_types.get_block(i);
    if (bl)
    {
      std::string name = ge_gen_types.get_name(i);
      std::string wood_tex_name = bl->get_string("wood_tex_name","wood");
      std::string leaf_tex_name = bl->get_string("leaf_tex_name","leaf");
      GETreeParameters *params = new GETreeParameters();
      params->load_from_blk(*bl);
      TreeTypeData type = TreeTypeData(id,params,wood_tex_name,leaf_tex_name);
      type.generator_name = "ge_gen";
      tree_types.emplace(name,type);
      id++;
    }
  }

  for (int i=0;i<my_gen_types.size();i++)
  {
    Block *bl = my_gen_types.get_block(i);
    if (bl)
    {
      std::string name = my_gen_types.get_name(i);
      std::string wood_tex_name = bl->get_string("wood_tex_name","wood");
      std::string leaf_tex_name = bl->get_string("leaf_tex_name","leaf");
      TreeStructureParameters *params = new TreeStructureParameters();
      params->load_from_blk(*bl);
      TreeTypeData type = TreeTypeData(id,params,wood_tex_name,leaf_tex_name);
      type.generator_name = "my_gen";
      tree_types.emplace(name,type);
      id++;
    }
  }
}
void base_init()
{
  glewInit();
  Tiny::view.lineWidth = 1.0f;
  Tiny::window("Procedural Tree", appContext.WIDTH, appContext.HEIGHT);
  Tiny::event.handler = [&]()
  { eventHandler(appContext, Tiny::event); };
  textureManager = TextureManager("./resources/textures/");
}
void init_render(WorldRenderer &worldRenderer)
{
  Block render_settings;
  worldRenderer.init(appContext.WIDTH, appContext.HEIGHT, render_settings);
}

void prepare_global_ggd_from_settings(SceneGenerator::SceneGenerationContext &ctx)
{
  for (auto p : ctx.tree_types)
  {
    ctx.global_ggd.types.push_back(p.second);
  }
}

int parser_main(int argc, char *argv[])
{
    base_init();
    Scene scene;
    SceneGenerator::SceneGenerationContext sceneGenerationContext;
    BlkManager man;
    man.load_block_from_file("scene_generation_settings.blk", sceneGenerationContext.settings);
    sceneGenerationContext.scene = &scene;
    load_tree_types(sceneGenerationContext.tree_types);
    prepare_global_ggd_from_settings(sceneGenerationContext);
    SceneGenerator sceneGen = SceneGenerator(sceneGenerationContext);
    sceneGen.create_scene_auto();

    bool need_render = true;
    if (need_render)
    {
        WorldRenderer worldRenderer;
        Block render_settings;
        worldRenderer.init(appContext.WIDTH, appContext.HEIGHT, render_settings);
        worldRenderer.set_heightmap(*scene.heightmap);
        worldRenderer.set_grove(scene.grove, sceneGenerationContext.global_ggd);
        Tiny::view.pipeline = [&]()
        {
            worldRenderer.set_forced_LOD(appContext.forced_LOD);
            worldRenderer.render(1,appContext.camera);
        };
        Tiny::view.interface = [&]()
        { 
        };

        Tiny::loop([&]() {});
        {
        };
    
        Tiny::quit();
    }

    return 0;
}

}