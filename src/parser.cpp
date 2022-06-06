#include "parser.h"
#define GLEW_EXPERIMENTAL
#include "tinyEngine/TinyEngine.h"
#include "glm/trigonometric.hpp"
#include "tinyEngine/image.h"
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

  void print_size(GrovePacked &grove)
  {
    int joints = 0;
    int leaves = 0;
    int matrixes = 0;
    for (auto &br : grove.instancedBranches)
    {
      matrixes += br.IDA.transforms.size();
      for (auto &bid : br.branches)
      {
        PackedBranch &b = grove.instancedCatalogue.get(bid);
        joints += b.joints.size();
        leaves += b.leaves.size();
      }
    }
    int verts = 4*joints + 4*leaves;
    int polys = 4*joints + 2*leaves;
    int bytes_per_vert = sizeof(float)*(3 + 3 + 2);
    int model_size = verts*bytes_per_vert + polys * 3 *sizeof(int);
    int mat_size = matrixes*sizeof(int)*(1+1+3+3+16);
    debug("Scene size:\n");
    debug("%.1fk joints\n", joints*1e-3);
    debug("%.1fk leaves\n", leaves*1e-3);
    debug("%.1fk vertices\n", verts*1e-3);
    debug("%.1fk polygons\n", polys*1e-3);
    debug("%d instance matrices\n", matrixes);
    debug("%.2f Mb models\n", model_size*1e-6);
    debug("%.2f Mb matrices\n", mat_size*1e-6);
    debug("%.2f Mb total\n", (model_size+ mat_size)*1e-6);
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
      sceneGen.set_biome_round(glm::vec2(0,0),gen_settings.get_double("forest_size",200),gen_settings.get_string("biome", "mixed_forest"));
      Block *objects = gen_settings.get_block("objects");
      if (objects)
      {
        for (int i = 0; i < objects->size(); i++)
        {
          Block *b = objects->get_block(i);
          if (b)
            sceneGen.add_object_blk(*b);
        }        
      }
      sceneGenerationContext.objects_bvh.rebuild();
      sceneGen.create_scene_auto();
      /*
      //sceneGen.set_biome_round(glm::vec2(0,0),100,"bush");
      int sz = ceil(sqrt((float)demo_mode_trees_cnt));
      int cnt = 0;
      float dist = 75;
      int t_id = metainfoManager.get_tree_type_id_by_name("medium_oak_simplified");
      int b_id = metainfoManager.get_tree_type_id_by_name("ref30.jpg_selected_params_9");
      int l_id = metainfoManager.get_tree_type_id_by_name("baobab");

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


      //sceneGen.set_biome_round(glm::vec2(0,0),500,"meadow");
      //sceneGen.set_biome_round(glm::vec2(0,0),gen_settings.get_double("forest_size",200),gen_settings.get_string("biome", "mixed_forest"));
      //sceneGen.set_biome_round(glm::vec2(0,0),gen_settings.get_double("forest_size",200),gen_settings.get_string("biome", "pine_forest"));
      //sceneGen.set_biome_round(glm::vec2(0,100),50,"bush");
      //sceneGenerationContext.biome_map.save_as_image();
      //
      Block objs;
      int p_cnt = 1;
      for (int i=0;i<p_cnt;i++)
      {
        for (int j=0;j<2;j++)
        {
          glm::vec pos = glm::vec2(150 * (i + urand(-0.2, 0.2)), 150 * (-1 + 2*j));
          sceneGen.plant_tree(pos, l_id);
        }
      }
      for (int i = 0; i < 0; i++)
      {
        Block *chb = new Block();
        chb->add_string("name", "farm_1");
        chb->add_bool("on_terrain", true);
        glm::vec3 pos = glm::vec3(100 * (i + urand(0, 0.5)), 0, 50);
        //sceneGen.plant_tree(glm::vec2(pos.x, 56 + pos.z), t_id);
        //sceneGen.plant_tree(glm::vec2(67 + pos.x, 10 + pos.z), l_id);
        glm::vec3 size = glm::vec3(125);
        pos.y = sceneGenerationContext.scene->heightmap->get_height(pos) - 0.4 * size.y;
        chb->add_mat4("transform", glm::scale(
                                       glm::rotate(
                                           glm::translate(glm::mat4(1.0f), pos),
                                           0*(float)urand() * PI, glm::vec3(0, 1, 0)),
                                       size));
        objs.add_block("obj", chb);
      }
      cnt = 0;
      for (int i = 0; i < cnt * cnt; i++)
      {
        Block *chb = new Block();
        chb->add_string("name", "stone_1");
        chb->add_bool("on_terrain", false);
        glm::vec3 pos = glm::vec3(100 * (i / cnt - cnt / 2 + urand(0, 0.67)), 0, 100 * (i % cnt - cnt / 2 + urand(0, 0.67)));
        glm::vec3 size = glm::vec3(urand(15, 40));
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
      man.save_block_to_file("objs.blk", objs);
      
      //print_size(scene.grove);
    */
    }
    else if (sandbox)
    {
      sandbox_main(argc, argv, &scene);
    }
    if (save_to_hydra)
    {
      HydraSceneExporter hExp;
      Block export_settings;
      glm::vec3 camera_pos = glm::vec3(-300, 100, -300);
      glm::vec3 camera_dir = glm::vec3(1, 0, 0.15);
      export_settings.add_vec3("camera_look_at", glm::vec3(0,130,0));
      export_settings.add_vec3("camera_pos", camera_pos);
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
        worldRenderer.set_render_mode(appContext.render_mode);
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