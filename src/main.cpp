
#define GLEW_EXPERIMENTAL
#include "tinyEngine/TinyEngine.h"
#include "generated_tree.h"
#include "glm/trigonometric.hpp"
#include "tinyEngine/helpers/image.h"
#include "tinyEngine/helpers/color.h"
#include "tinyEngine/helpers/helper.h"
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include <string>
#include "tinyEngine/camera.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "tinyEngine/utility.h"
#include "grove.h"
#include "tinyEngine/save_utils/config.h"
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include "terrain.h"
#include "tinyEngine/shadow.h"
#include "grove_packer.h"
#include "proctree.h"
#include "app.h"
#include "grass_renderer.h"
#include "tinyEngine/utility/deffered_target.h"
#include "tinyEngine/ambient_occlusion.h"
#include "tinyEngine/utility/cubemap.h"
#include "tinyEngine/utility/python_interaction.h"
#include "tinyEngine/gltf_utils/general_gltf_writer.h"
#include <thread>
#include "parameter_selection.h"
#include "tinyEngine/save_utils/blk.h"
#include "clustering/clustering_benchmark.h"
#include "simple_generator.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
AppContext appContext;
ShadowMap shadowMap;
DefferedTarget defferedTarget;
GrovePacker packer;
AbstractTreeGenerator *gen = nullptr;
DebugVisualizer *debugVisualizer = nullptr;
GrovePacked grove;
GroveGenerationData ggd;
GroveRenderer *GR = nullptr;
ClusteringBenchmark benchmark;
struct RenderData
{
  bool regenerate_shadows = true;

  HBAORenderer *hbaoRenderer;
  Cubemap *cubemap;
  PostFx *defferedLight;
  PostFx *startScreenShader;
  Shader *defaultShader;
  Shader *debugShader;

  Heightmap *heightmap;
  GroveRenderer *groveRenderer;
  HeightmapTex *heightmapTex;
  GrassRenderer *grassRenderer;
  TerrainRenderer *terrainRenderer;

  void clear();
} data;

bool generation_needed = false;
bool saving_needed = false;
bool loading_needed = false;
bool print_perf = false;
bool only_gen = false;
bool visualize_voxels = false;
bool visualize_clusters = false;
bool visualize_initial_voxels = false;
bool statistics_run = false;
bool gltf_export = false;
bool need_initialization = true;
bool parameter_selection = false;
bool clustering_benchmark = false;
bool prepare_dataset = false;
std::string generator_name = "default";
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
std::string save_dataset_path = ".";
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
      visualize_clusters = true;
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
        logerr("possible names : default, proctree, simple");
      }
      k += 3;
    }
    else if (std::string(argv[k]) == "-prepare_dataset")
    {
      prepare_dataset = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -prepare_dataset");
        return 1;
      }
      else
        save_dataset_path = argv[k + 1];
      k += 2;
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
bool prepare_dictory(std::string &save_path)
{
  bool status = true;
  try
  {
    if (boost::filesystem::exists(save_path))
    {
      if (boost::filesystem::is_directory(save_path))
      {
        printf("replacing previous save\n");
        boost::filesystem::remove_all(save_path);
      }
      else
      {
        logerr("path %s represents existing file. Can not save here", save_path.c_str());
        status = false;
      }
    }
    if (status)
    {
      boost::filesystem::create_directory(save_path);
      boost::filesystem::permissions(save_path, boost::filesystem::perms::all_all);
    }
  }
  catch (const std::exception &e)
  {
    status = false;
    std::cerr << e.what() << '\n';
  }

  return status;
}
void clear_current_grove()
{
  grove = GrovePacked();
  if (GR)
  {
    delete GR;
    GR = nullptr;
    data.groveRenderer = nullptr;
  }
  packer = GrovePacker();
}
void generate_grove()
{
  ::Tree *trees = new ::Tree[ggd.trees_count];
  gen->create_grove(ggd, trees, *data.heightmap);
  logerr("%d branches",trees[0].branchHeaps[1]->branches.size());
  if (prepare_dataset)
  {
    bool status = prepare_dictory(save_dataset_path);
    if (status)
      packer.add_trees_to_grove_prepare_dataset(ggd, grove, trees, data.heightmap, save_dataset_path);
    else
    {
      logerr("unable to create directory to save dataset. Exiting.");
    }
  }
  else
    packer.add_trees_to_grove(ggd, grove, trees, data.heightmap, visualize_clusters);
  delete[] trees;
}
void generate_single_tree(ParametersSet *par, GrovePacked &res)
{
    GrovePacker packer;
    GroveGenerationData tree_ggd;
    DistributionData dd;
    distibutionGenerator.d = &dd;

    tree_ggd.trees_count = 1;
    TreeTypeData type = ggd.types[0];
    type.params = par;
    tree_ggd.types = {type};
    tree_ggd.synts_count = 0;
    tree_ggd.synts_precision = 1;
    tree_ggd.pos = glm::vec3(0,0,0);
    tree_ggd.size = glm::vec3(150,150,250);
    tree_ggd.obstacles = {};
    tree_ggd.clustering_max_individual_distance = 0.0;
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    Tree single_tree;
    gen->create_grove(tree_ggd, &single_tree, *data.heightmap);
    packer.add_trees_to_grove(ggd, res, &single_tree, data.heightmap, false);
    print_alloc_info();
    distibutionGenerator.d = nullptr;
    dd.clear();
}
void generate_grove_renderer(GrovePacked *source, GroveGenerationData *gen_data)
{
  std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
  if (pres == GroveRenderer::Precision::LOW)
    LODs_dists.back() = -10;
  if (data.groveRenderer)
  {
    delete data.groveRenderer;
  }
  data.groveRenderer = new GroveRenderer(source, gen_data, 5, LODs_dists, print_perf, pres);
  GR = data.groveRenderer;
}
int base_initialization()
{
  //base initialization
  glewInit();
  Tiny::view.lineWidth = 1.0f;
  Tiny::window("Procedural Tree", appContext.WIDTH, appContext.HEIGHT);
  Tiny::event.handler = [&]()
  { eventHandler(appContext, Tiny::event); };
  textureManager = TextureManager("./resources/textures/");
  data.startScreenShader = new PostFx("simple_render.fs");

  //PythonHelper ph;
  //ph.test();
  return -1;
}
int full_initialization()
{
  config.load_config();
  config.load_ggds();

  appContext.camera.pos = glm::vec3(-300, 70, 0);
  appContext.light.dir = glm::normalize(glm::vec3(-0.2, 0.5, -0));
  appContext.light.color = glm::vec3(0.99, 0.9, 0.7);
  appContext.light.intensity = 1;
  appContext.light.ambient_q = 0.37;
  appContext.light.diffuse_q = 0.63;
  appContext.light.specular_q = 0.0;
  appContext.light.has_shadow_map = true;
  appContext.light.shadow_map_size = glm::vec2(4096, 4096);
  shadowMap.create(4096, 4096);
  defferedTarget.create(2 * 1920, 2 * 1080);
  defferedTarget.set_clear_color(glm::vec4(0.0, 0.0, 0.0, 0.0));

  data.hbaoRenderer = new HBAORenderer();
  data.hbaoRenderer->create(800, 600);
  data.cubemap = new Cubemap(1920, 1080);
  data.defferedLight = new PostFx("deffered_light.fs");
  data.debugShader = new Shader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  data.defaultShader = new Shader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), data.defaultShader);

  srand(1);

  data.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(2000, 2000), 5);
  if (generator_name == "proctree")
    gen = new Proctree::ProctreeGenerator();
  else if (generator_name == "simple")
    gen = new SimpleTreeGenerator();
  else
    gen = new mygen::TreeGenerator();
  data.heightmap->random_generate(0, 0, 0);

  if (generation_needed)
  {
    ggd = config.get_ggd(grove_type_name);
    if (generator_name == "proctree")
    {
      ggd.types[0].params = new Proctree::Properties();
    }
    if (parameter_selection)
    {
      std::function<void(ParametersSet *, GrovePacked &)> _generate = generate_single_tree;
      ParameterSelector sel(_generate);
      ParametersSet *start = ggd.types[0].params;
      Quality imp_qual = ggd.impostor_quality;
      ggd.impostor_quality = Quality::ULTRALOW;
      //ggd.bill_1_quality = Quality::ULTRALOW;
      //ggd.bill_2_quality = Quality::ULTRALOW;
      sel.select(start,parameter_selector_name);
      ggd.impostor_quality = imp_qual;
    }
    
    if (statistics_run)
    {
      ggd.clustering_max_individual_distance = statRunLaunchParams.max_ind_dist;
      ggd.trees_count = statRunLaunchParams.trees;
    }
    if (clustering_benchmark)
    {
      benchmark.perform_benchmark(clustering_benchmark_path, gen, ggd, data.heightmap);
      if (benchmark.grove_count() > 0)
      {
        appContext.benchmark_grove_needed = 0;
      }
    }
    else
    {
      generate_grove();
    }
    if (visualize_initial_voxels)
    {
      mygen::TreeGenerator *mygen_gen = dynamic_cast<mygen::TreeGenerator *>(gen);
      if (mygen_gen)
        debugVisualizer->visualize_light_voxels(mygen_gen->voxels);
    }
    /*LightVoxelsCube *voxx = new LightVoxelsCube(glm::vec3(0,0,0), glm::ivec3(8,8,8),1,3,2);
    for (int i=-8;i<=8;i++)
    for (int j=-8;j<=8;j++)
    for (int k=-8;k<=8;k++)
    voxx->set_occluder_trilinear(glm::vec3(i,j,k),abs(10*i));
    voxx->prepare_mip_levels();
    for (int i=0;i<0;i++)
    {
      logerr("visualize voxels");
      if (test_voxels_cube[0])
      {
        auto *voxels = test_voxels_cube[0];
        debugVisualizer->visualize_light_voxels(voxels,
                            voxels->get_center() - voxels->get_voxel_size()*glm::vec3(voxels->get_vox_sizes()),
                            voxels->get_voxel_size()*(2.0f*glm::vec3(voxels->get_vox_sizes()) + glm::vec3(1)),
                            powf(2,i)*glm::vec3(voxels->get_voxel_size()),
                            powf(2,i)*0.85f*voxels->get_voxel_size(),
                            0.01,glm::vec3(0,100*(i+1),0),glm::vec3(1,1,1),i);
      }
    }*/
  }
  if (saving_needed)
  {
    bool status = prepare_dictory(save_path);
    if (status)
    {
      std::string f_path = save_path + std::string("/grove.dat");
      FILE *f = fopen(f_path.c_str(), "wb");
      saver::set_textures_path(save_path);
      saver::save(f, grove);
      fclose(f);
    }
    else
    {
      logerr("error occured while saving to path %s", save_path.c_str());
    }
  }
  if (loading_needed)
  {
    try
    {
      if (boost::filesystem::exists(load_path))
      {
        if (boost::filesystem::is_directory(load_path))
        {
          grove = GrovePacked();
          struct stat sb;
          if (stat(load_path.c_str(), &sb) != 0)
          {
            mkdir(load_path.c_str(), 0777);
          }
          if (!S_ISDIR(sb.st_mode))
          {
            logerr("given load path \"%s\" is not a directory", load_path.c_str());
          }
          std::string f_path = load_path + std::string("/grove.dat");
          FILE *f = fopen(f_path.c_str(), "rb");
          saver::set_textures_path(load_path);
          saver::load(f, grove);
          ggd = config.get_ggd(grove.ggd_name);
          logerr("grove %s", grove.ggd_name.c_str());
          fclose(f);
        }
        else
        {
          logerr("given load path %s is not a directory", load_path.c_str());
        }
      }
      else
      {
        logerr("given load path %s does not exist", load_path.c_str());
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
    }
  }
  generate_grove_renderer(&grove, &ggd);

  for (int i = 0; i < ggd.obstacles.size(); i++)
    debugVisualizer->add_bodies(ggd.obstacles[i], 1);
  data.terrainRenderer = new TerrainRenderer(*data.heightmap, glm::vec3(0, 0, 0), glm::vec2(2500, 2500), glm::vec2(25, 25));

  data.heightmapTex = new HeightmapTex(*data.heightmap, 2048, 2048);
  data.grassRenderer = new GrassRenderer();

  if (gltf_export)
  {
    gltf::GeneralGltfWriter ggw;
    ggw.add_packed_grove(grove, ggd);
    ggw.add_model(data.terrainRenderer->flat_terrain);
    ggw.convert_to_gltf("scene");
    ggw.clear();
    debug("Scene exported to .gltf format\n");
  }
  if (only_gen)
  {
    debug("Grove successfully generated. Exiting.");
    return 0;
  }

  appContext.renderMode = RenderMode::Rendering;
  Texture tex = textureManager.get("reference_tree_test");
  return -1;
}

void simple_render_pipeline()
{
  auto &ctx = appContext;

  ctx.fpsCounter.tick();
  if (ctx.fpsCounter.get_frame_n() % 100 == 0 && print_perf)
  {
    fprintf(stderr, "FPS: %4.1f\n", ctx.fpsCounter.get_average_fps());
  }
  if (data.regenerate_shadows)
  {
    //shadows pass
    data.regenerate_shadows = false;
    shadowMap.use(ctx.light);
    glm::mat4 sh_viewproj = shadowMap.get_transform();
    if (data.groveRenderer)
    {
      data.groveRenderer->render(data.groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), ctx.camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                ctx.light, ctx.groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (data.terrainRenderer)
    {
      data.terrainRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(),
                                  0, ctx.camera.pos, ctx.light, true);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                ctx.camera.pos, *data.heightmapTex, ctx.light, true);
    }
    shadowMap.start_trans_pass();
    if (data.groveRenderer)
    {
      data.groveRenderer->render(data.groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), ctx.camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                ctx.light, ctx.groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                ctx.camera.pos, *data.heightmapTex, ctx.light, true);
    }
    shadowMap.finish_trans_pass();

    shadowMap.blur();
  }
  //Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  defferedTarget.target();

  if (ctx.render_mode <= ctx.DEBUG_RENDER_MODE)
  {
    data.debugShader->use();
    data.debugShader->uniform("projectionCamera", ctx.projection * ctx.camera.camera());
    if (ctx.render_mode == ctx.DEBUG_RENDER_MODE)
    {
      data.debugShader->texture("tex", textureManager.get(ctx.debug_tex));
      data.debugShader->texture("tex_arr", textureManager.get(ctx.debug_tex));
      data.debugShader->uniform("need_tex", true);
      data.debugShader->uniform("need_arr_tex", false);
      data.debugShader->uniform("need_coord", false);
      data.debugShader->uniform("slice", 0);
    }
    else if (ctx.render_mode == ctx.ARRAY_TEX_DEBUG_RENDER_MODE)
    {
      data.debugShader->texture("tex", textureManager.get_arr(ctx.debug_tex));
      data.debugShader->texture("tex_arr", textureManager.get_arr(ctx.debug_tex));
      data.debugShader->uniform("need_tex", false);
      data.debugShader->uniform("need_arr_tex", true);
      data.debugShader->uniform("need_coord", false);
      data.debugShader->uniform("slice", ctx.debug_layer);
    }
  }
  else
  {
    //depth prepass
    /*tr.render(ctx.projection * ctx.camera.camera(),shadowMap.get_transform(),shadowMap.getTex(),
                ctx.camera.pos,ctx.light, true);

      glClearColor(clearcolor.x, clearcolor.y, clearcolor.z, 1.0f);*/
    //color pass
    if (data.terrainRenderer)
    {
      data.terrainRenderer->render(ctx.projection, ctx.camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                  ctx.camera.pos, ctx.light);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(ctx.projection, ctx.camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                ctx.camera.pos, *data.heightmapTex, ctx.light);
    }
    if (ctx.render_mode != 2 && data.groveRenderer)
    {
      data.groveRenderer->render(ctx.forced_LOD, ctx.projection, ctx.camera.camera(), ctx.camera,
                                 glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT), ctx.light,
                                 ctx.groveRendererDebugParams, shadowMap.get_transform(), 0 * shadowMap.getTex());
    }
    if (debugVisualizer)
      debugVisualizer->render(ctx.projection * ctx.camera.camera(), ctx.render_mode);
  }

  //postfx
  /*uniform vec3 dir_to_sun;
uniform vec3 camera_pos;
uniform vec3 ambient_diffuse_specular;
uniform vec3 light_color;
uniform vec2 sts_inv;
uniform sampler2D shadowMap;
uniform bool need_shadow;
uniform mat4 shadow_mat;*/
  //hbaoRenderer.render(ctx,defferedTarget.get_view_pos());

  data.cubemap->render(ctx.projection, ctx.camera.camera(), ctx.camera);
  glm::vec3 ads = glm::vec3(ctx.light.ambient_q, ctx.light.diffuse_q, ctx.light.specular_q);
  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  data.defferedLight->use();
  data.defferedLight->get_shader().texture("colorTex", defferedTarget.get_color());
  data.defferedLight->get_shader().texture("normalsTex", defferedTarget.get_normals());
  data.defferedLight->get_shader().texture("viewPosTex", defferedTarget.get_view_pos());
  data.defferedLight->get_shader().texture("worldPosTex", defferedTarget.get_world_pos());
  data.defferedLight->get_shader().texture("aoTex", data.hbaoRenderer->get_tex());
  data.defferedLight->get_shader().texture("shadowMap", shadowMap.getTex());
  data.defferedLight->get_shader().texture("cubeTex", data.cubemap->get_tex());
  data.defferedLight->get_shader().uniform("dir_to_sun", ctx.light.dir);
  data.defferedLight->get_shader().uniform("camera_pos", ctx.camera.pos);
  data.defferedLight->get_shader().uniform("ambient_diffuse_specular", ads);
  data.defferedLight->get_shader().uniform("light_color", ctx.light.color);
  data.defferedLight->get_shader().uniform("need_shadow", shadowMap.getTex() != 0);
  data.defferedLight->get_shader().uniform("shadow_mat", shadowMap.get_transform());
  data.defferedLight->get_shader().uniform("sts_inv", 1.0f / ctx.light.shadow_map_size);
  data.defferedLight->render();
}

void start_screen_pipeline()
{
  Texture start_screen = textureManager.get("start_screen");
  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  data.startScreenShader->use();
  data.startScreenShader->get_shader().texture("tex", start_screen);
  data.startScreenShader->get_shader().uniform("tex_transform", glm::vec4(0, 1, 1, -1));
  data.startScreenShader->render();
}

int main(int argc, char *argv[])
{
  int status = parse_arguments(argc, argv);
  if (status >= 0)
    return status;

  status = base_initialization();
  if (status >= 0)
    return status;

  if (need_initialization && appContext.renderMode != RenderMode::Rendering)
  {
    status = full_initialization();
    if (status >= 0)
      return status;
  }

  Tiny::view.pipeline = [&]()
  {
    if (need_initialization && appContext.renderMode != RenderMode::Rendering)
    {
      full_initialization();
      //std::thread t1([&]{full_initialization();});
      //t1.join();
    }
    if (appContext.regeneration_needed)
    {
      clear_current_grove();
      generate_grove();
      generate_grove_renderer(&grove, &ggd);
      //print_alloc_info();
      appContext.regeneration_needed = false;
    }
    else if (appContext.add_generation_needed)
    {
      //print_alloc_info();
      generate_grove();
      //print_alloc_info();
      generate_grove_renderer(&grove, &ggd);
      //print_alloc_info();
      appContext.add_generation_needed = false;
    }
    if (appContext.benchmark_grove_current != appContext.benchmark_grove_needed &&
        appContext.benchmark_grove_needed >= 0)
    {
      int g = appContext.benchmark_grove_needed % benchmark.grove_count();
      appContext.benchmark_grove_current = appContext.benchmark_grove_needed;
      generate_grove_renderer(&(benchmark.get_grove(g)), &ggd);
    }
    if (appContext.renderMode == RenderMode::StartingScreen)
      start_screen_pipeline();
    else if (appContext.renderMode == RenderMode::Rendering)
      simple_render_pipeline();
  };

  Tiny::view.interface = [&]()
  {
    bool p = true;
    //ImGui::ShowDemoWindow(&p);
  };

  Tiny::loop([&]() {});
  {
    for (int i = 0; i < ggd.obstacles.size(); i++)
    {
      delete ggd.obstacles[i];
      ggd.obstacles[i] = nullptr;
    }
    ggd.obstacles.clear();
  }

  data.clear();
  distibutionGenerator.clear();
  Tiny::quit();

  return 0;
}

void RenderData::clear()
{
  regenerate_shadows = true;
  if (hbaoRenderer)
    delete hbaoRenderer;
  if (cubemap)
    delete cubemap;
  if (defferedLight)
    delete defferedLight;
  if (startScreenShader)
    delete startScreenShader;
  if (defaultShader)
    delete defaultShader;
  if (debugShader)
    delete debugShader;
  if (groveRenderer)
    delete groveRenderer;
  if (heightmapTex)
    delete heightmapTex;
  if (grassRenderer)
    delete grassRenderer;
  if (terrainRenderer)
    delete terrainRenderer;
  if (heightmap)
    delete heightmap;
}
