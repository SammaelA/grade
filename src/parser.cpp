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
#include "generation/metainfo_manager.h"
#include "sandbox.h"

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
  bool sandbox = false;
  int demo_mode_trees_cnt = 0;
  int demo_mode_patch_size = 5;
  bool debug_bvh = false;
  bool no_init = false;
  int parse_arguments(int argc, char *argv[])
  {
    int k = 1;
    while (k < argc)
    {
      if (std::string(argv[k]) == "-settings")
      {
        if (argc > k + 2 && std::string(argv[k + 1]) == "=")
        {
          settings_block = std::string(argv[k + 2]);
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
      else if (std::string(argv[k]) == "-no_init")
      {
        no_init = true;
        k++;
      }
      else if (std::string(argv[k]) == "-hydra")
      {
        save_to_hydra = true;
        if (argc > k + 2 && std::string(argv[k + 1]) == "=")
        {
          hydra_scene_dir = std::string(argv[k + 2]);
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
          int n = std::stoi(argv[k + 1]);
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
          int n = std::stoi(argv[k + 1]);
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
      else if (std::string(argv[k]) == "-debug_bvh")
      {
        debug_bvh = true;
        k++;
      }
      else if (std::string(argv[k]) == "-no_debug")
      {
        debug_level = 1000;
        k++;
      }
      else if (std::string(argv[k]) == "-sandbox")
      {
        sandbox = true;
        k++;
      }
      else if (std::string(argv[k]) == "-debug")
      {
        if (argc > k + 1)
        {
          int n = std::stoi(argv[k + 1]);
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
      else if (std::string(argv[k]) == "-random_seed")
      {
        srand(time(nullptr));
        k++;
      }
      else
      {
        logerr("unknows argument %s", argv[k]);
        k++;
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
    man.load_block_from_file("resources.blk", textures_list);
    textureManager = TextureManager("./resources/textures/", textures_list);
    ModelLoader::load_default_blk();
  }
  void init_render(WorldRenderer &worldRenderer)
  {
    Block render_settings;
    worldRenderer.init(appContext.WIDTH, appContext.HEIGHT, render_settings);
  }

  void demo_scene_ctx(Block &settings)
  {
    demo_mode_patch_size = MIN(demo_mode_patch_size, demo_mode_trees_cnt);
    float density = MAX(demo_mode_patch_size, 2);
    int d_sqrt = sqrt(density);
    glm::vec2 cell_size(30 * (d_sqrt + 1), 30 * (d_sqrt + 1));
    int patches_cnt = ceil(demo_mode_trees_cnt / density);
    int patches_x = sqrt(patches_cnt);
    patches_x = MAX(patches_x, 1);
    int patches_y = patches_x;
    int max_trees_per_patch = 2 * ((float)demo_mode_trees_cnt / patches_cnt);
    settings.set_vec2("cell_size", cell_size);
    settings.set_vec2("scene_size", glm::vec2(cell_size.x * patches_x, cell_size.y * patches_y));
    settings.set_int("max_trees_per_patch", max_trees_per_patch);
    settings.set_int("fixed_patches_count", patches_cnt);
    settings.set_double("patches_density", 0);
  }

  int parser_main(int argc, char *argv[])
  {
    parse_arguments(argc, argv);
    if (sandbox && no_init)
    {
      sandbox_main(argc, argv, nullptr);
      exit(0);
    }
    base_init();
    Scene scene;
    SceneGenerator::SceneGenerationContext sceneGenerationContext;
    sceneGenerationContext.scene = &scene;
    BlkManager man;
    Block gen_settings;
    parse_arguments(argc, argv);
    man.load_block_from_file(settings_block, gen_settings);
    if (demo_mode)
    {
      demo_scene_ctx(gen_settings);

      SceneGenerator sceneGen = SceneGenerator(sceneGenerationContext);

      sceneGen.init_scene(gen_settings);
      sceneGen.create_heightmap_simple_auto();
      sceneGen.set_default_biome("meadow");
      //sceneGen.set_biome_round(glm::vec2(0,0),100,"bush");
      int sz = ceil(sqrt((float)demo_mode_trees_cnt));
      int cnt = 0;
      float dist = 75;
      int t_id = metainfoManager.get_tree_type_id_by_name("sphere_tree");
      int b_id = metainfoManager.get_tree_type_id_by_name("apple");
      int l_id = metainfoManager.get_tree_type_id_by_name("large_oak");

      int ids[3] = {t_id, b_id, l_id};
      for (int i = 0; i < sz; i++)
      {
        for (int j = 0; j < sz; j++)
        {
          if (cnt < demo_mode_trees_cnt)
          {
            float tx = dist * (j - sz / 2 + urand());
            float ty = dist * (i - sz / 2 + urand());
            sceneGen.plant_tree(glm::vec2(tx, ty), ids[cnt % 3]);
            cnt++;
          }
        }
      }

      //sceneGen.set_biome_round(glm::vec2(0,0),50,"mixed_forest");
      //sceneGen.set_biome_round(glm::vec2(0,100),50,"bush");
      sceneGenerationContext.biome_map.save_as_image();
      //
      Block objs;
      cnt = 0;
      for (int i = 0; i < cnt * cnt; i++)
      {
        Block *chb = new Block();
        chb->add_string("name", "stone_1");
        chb->add_bool("on_terrain", false);
        glm::vec3 pos = glm::vec3(50 * (i / cnt - cnt / 2 + urand(-1, 1)), 0, 50 * (i % cnt - cnt / 2 + urand(-1, 1)));
        glm::vec3 size = glm::vec3(urand(10, 30));
        pos.y = sceneGenerationContext.scene->heightmap->get_height(pos) - 0.4 * size.y;
        chb->add_mat4("transform", glm::scale(
                                       glm::rotate(
                                           glm::translate(glm::mat4(1.0f), pos),
                                           (float)urand() * PI, glm::vec3(0, 1, 0)),
                                       size));
        objs.add_block("obj", chb);
      }

      for (int i = 0; i < objs.size(); i++)
      {
        Block *b = objs.get_block(i);
        if (b)
          sceneGen.add_object_blk(*b);
      }
      //
      sceneGenerationContext.objects_bvh.rebuild();
      sceneGen.create_scene_auto();
    }
    else if (sandbox)
    {
      sandbox_main(argc, argv, &scene);
    }
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
      GroveGenerationData ggd;
      ggd.types = metainfoManager.get_all_tree_types();
      worldRenderer.init(appContext.WIDTH, appContext.HEIGHT, render_settings);
      worldRenderer.set_heightmap(*scene.heightmap);
      worldRenderer.set_grass(scene.grass);
      worldRenderer.set_grove(scene.grove, ggd);
      auto &voxels = debugTransferData.debug_voxels;
      for (auto *vox : voxels)
      {
        if (vox)
          worldRenderer.set_voxels_debug(*vox);
      }
      worldRenderer.add_instanced_models(scene.instanced_models);

      if (debug_bvh)
      {
        auto func = [&](const std::pair<AABB, uint64_t> &p)
        {
          worldRenderer.add_aabb_debug(p.first);
        };
        sceneGenerationContext.objects_bvh.iterate_over_intersected_bboxes(AABB(glm::vec3(-1e9, -1e9, -1e9),
                                                                                glm::vec3(1e9, 1e9, 1e9)),
                                                                           func, false);
      }
      Tiny::view.pipeline = [&]()
      {
        worldRenderer.set_resolution(Tiny::view.WIDTH, Tiny::view.HEIGHT);
        worldRenderer.set_forced_LOD(appContext.forced_LOD);
        worldRenderer.render(1, appContext.camera);
      };
      Tiny::view.interface = [&]() {
      };

      Tiny::loop([&]() {});
      {};

      Tiny::quit();
    }

    return 0;
  }

}