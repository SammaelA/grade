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
#include "hydra_utils/hydra_scene_exporter.h"
#include "graphics_utils/debug_transfer.h"

namespace parser
{
//View Tiny::view;   //Window and Interface  (Requires Initialization)
//Event Tiny::event; //Event Handler
//Audio Tiny::audio; //Audio Processor       (Requires Initialization)
AppContext appContext;

bool render_needed = false;
bool save_to_hydra = false;
std::string settings_block = "scene_generation_settings.blk";
std::string hydra_scene_dir = "vegetation_scene";
bool demo_mode = false;
int demo_mode_trees_cnt = 0;
int demo_mode_patch_size = 5;
int parse_arguments(int argc, char *argv[])
{
  int k = 1;
  while (k < argc)
  {
    if (std::string(argv[k]) == "-settings")
    {
      if (argc > k + 2 && std::string(argv[k+1]) == "=")
      {
        settings_block = std::string(argv[k+2]);
        k += 3;
      }
      else
      {
        logerr("use \"settings = <settings_file.blk>\"");
        k++;
      }
    }
    else if (std::string(argv[k]) == "-render")
    {
      render_needed = true;
      k++;
    }
    else if (std::string(argv[k]) == "-no_render")
    {
      render_needed = false;
      k++;
    }
    else if (std::string(argv[k]) == "-hydra")
    {
      save_to_hydra = true;
      if (argc > k + 2 && std::string(argv[k+1]) == "=")
      {
        hydra_scene_dir = std::string(argv[k+2]);
        k += 3;
      }
      else
      {
        k++;
      }
    }
    else if (std::string(argv[k]) == "-demo")
    {
      if (argc > k + 1)
      {
        int n = std::stoi(argv[k+1]);
        if (n <= 0)
        {
          logerr("use \"demo <trees_count>\"");
        }
        else
        {
          demo_mode = true;
          demo_mode_trees_cnt = n;
        }
        k += 2;
      }
      else
      {
        logerr("use \"demo <trees_count>\"");
        k++;
      }
    }
    else if (std::string(argv[k]) == "-patch_size")
    {
      if (argc > k + 1)
      {
        int n = std::stoi(argv[k+1]);
        if (n <= 0)
        {
          logerr("use \"patch_size <trees_count_in_patch>\"");
        }
        else
        {
          demo_mode_patch_size = n;
        }
        k += 2;
      }
      else
      {
        logerr("use \"patch_size <trees_count_in_patch>\"");
        k++;
      }
    }
    else if (std::string(argv[k]) == "-debug_small_voxels")
    {
      debugTransferSettings.save_small_voxels_count = 100000;
      k++;
    }
    else if (std::string(argv[k]) == "-debug_original_voxels")
    {
      debugTransferSettings.save_detailed_voxels_count = 100000;
      k++;
    }
    else if (std::string(argv[k]) == "-no_debug")
    {
      debug_level = 1000;
      k++;
    }
    else if (std::string(argv[k]) == "-debug")
    {
      if (argc > k + 1)
      {
        int n = std::stoi(argv[k+1]);
        if (n <= 0)
        {
          logerr("use \"debug <debug_type>\"");
        }
        else
        {
          debug_level = n;
        }
        k += 2;
      }
      else
      {
        logerr("use \"debug <debug_type>\"");
        k++;
      }
    }
    else
    {
      logerr("unknows argument %s",argv[k]);
      k++;
    }
  }
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

void load_grass_types(std::map<std::string, GrassType> &grass_types)
{
  BlkManager man;
  Block grass_types_blk;
  int id = 0;
  man.load_block_from_file("grass_presets.blk",grass_types_blk);
  for (int i=0;i<grass_types_blk.size();i++)
  {
    Block *bl = grass_types_blk.get_block(i);
    std::string name = grass_types_blk.get_name(i);
    if (bl)
    {
      GrassType gt;
      gt.id = id;
      gt.load_from_blk(*bl);
      grass_types.emplace(name,gt);
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
  BlkManager man;
  Block textures_list;
  man.load_block_from_file("resources.blk",textures_list);
  textureManager = TextureManager("./resources/textures/", textures_list);
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

  ctx.global_ggd.obstacles.push_back(new Box(glm::vec3(-200,0,-20), glm::vec3(30,0,0), glm::vec3(0,0,30), glm::vec3(0,100,0)));
  ctx.global_ggd.obstacles.push_back(new Box(glm::vec3(-175,0,-50), glm::vec3(30,0,0), glm::vec3(0,0,30), glm::vec3(0,130,0)));
  ctx.global_ggd.obstacles.push_back(new Box(glm::vec3(-200,0,75), glm::vec3(30,0,0), glm::vec3(0,0,30), glm::vec3(0,90,0)));
  ctx.global_ggd.obstacles.push_back(new Box(glm::vec3(-230,0,-75), glm::vec3(30,0,0), glm::vec3(0,0,30), glm::vec3(0,110,0)));
  ctx.global_ggd.obstacles.push_back(new Box(glm::vec3(-270,0,50), glm::vec3(30,0,0), glm::vec3(0,0,30), glm::vec3(0,120,0)));
}

void demo_scene_ctx(SceneGenerator::SceneGenerationContext &sceneGenerationContext)
{
  demo_mode_patch_size = MIN(demo_mode_patch_size, demo_mode_trees_cnt);
  float density = MAX(demo_mode_patch_size, 2);
  int d_sqrt = sqrt(density);
  glm::vec2 cell_size(30*(d_sqrt + 1),30*(d_sqrt + 1));
  int patches_cnt = ceil(demo_mode_trees_cnt/density);
  int patches_x = sqrt(patches_cnt);
  patches_x = MAX(patches_x,1);
  int patches_y = patches_x;
  int max_trees_per_patch = 2*((float)demo_mode_trees_cnt/patches_cnt);
  sceneGenerationContext.settings.set_vec2("cell_size",cell_size);
  sceneGenerationContext.settings.set_vec2("scene_size",glm::vec2(cell_size.x*patches_x,cell_size.y*patches_y));
  sceneGenerationContext.settings.set_int("max_trees_per_patch",max_trees_per_patch);
  sceneGenerationContext.settings.set_int("fixed_patches_count", patches_cnt);
  sceneGenerationContext.settings.set_double("patches_density", 0);
}

int parser_main(int argc, char *argv[])
{
    base_init();
    Scene scene;
    SceneGenerator::SceneGenerationContext sceneGenerationContext;
    BlkManager man;

    parse_arguments(argc,argv);
    man.load_block_from_file(settings_block, sceneGenerationContext.settings);
    sceneGenerationContext.scene = &scene;
    load_tree_types(sceneGenerationContext.tree_types);
    load_grass_types(sceneGenerationContext.grass_types);
    prepare_global_ggd_from_settings(sceneGenerationContext);
    scene.tree_types = sceneGenerationContext.global_ggd.types;
    if (demo_mode)
      demo_scene_ctx(sceneGenerationContext);
    
    SceneGenerator sceneGen = SceneGenerator(sceneGenerationContext);
    sceneGen.create_scene_auto();

    if (save_to_hydra)
    {
      HydraSceneExporter hExp;
      Block export_settings;
      hExp.export_scene(hydra_scene_dir, scene, export_settings);
    }

    if (render_needed)
    {
        WorldRenderer worldRenderer;
        Block render_settings;
        worldRenderer.init(appContext.WIDTH, appContext.HEIGHT, render_settings);
        worldRenderer.set_heightmap(*scene.heightmap);
        worldRenderer.set_grass(scene.grass);
        worldRenderer.set_grove(scene.grove, sceneGenerationContext.global_ggd);
        for (auto *b : sceneGenerationContext.global_ggd.obstacles)
          worldRenderer.add_body_debug(b);
        auto &voxels = debugTransferData.debug_voxels;
        for (auto *vox : voxels)
        {
          if (vox)
            worldRenderer.set_voxels_debug(*vox);
        }
        Tiny::view.pipeline = [&]()
        {
          worldRenderer.set_resolution(Tiny::view.WIDTH, Tiny::view.HEIGHT);
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