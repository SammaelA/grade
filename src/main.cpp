
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
#include "tinyEngine/gltf_utils/general_gltf_writer.h"
#include <thread>
#include "parameter_selection.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
AppContext appContext;
ShadowMap shadowMap;
DefferedTarget defferedTarget;

Tree t[MAX_TREES];
mygen::TreeGenerator gen;
DebugVisualizer *debugVisualizer = nullptr;
GrovePacked grove;
GroveGenerationData ggd;
GroveRenderer *GR = nullptr;

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
bool visualize_initial_voxels = false;
bool statistics_run = false;
bool gltf_export = false;
bool need_initialization = true;
bool parameter_selection = false;
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
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_voxels")
    {
      visualize_voxels = true;
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
void clear_current_grove()
{
  grove = GrovePacked();
  if (GR)
  {
    delete GR;
    GR = nullptr;
    data.groveRenderer = nullptr;
  }
  gen.reset();
  for (int i=0;i<MAX_TREES;i++)
  {
    t[i].clear();
  }
}
void generate_grove()
{
    GrovePacker packer;
    gen.create_grove(ggd, t, *data.heightmap);
    logerr("%d branches",t[0].branchHeaps[1]->branches.size());
    packer.pack_grove(ggd, grove, *debugVisualizer, t, data.heightmap, visualize_voxels);
}
void generate_single_tree(TreeStructureParameters &par, GrovePacked &res)
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
    tree_ggd.clustering_max_individual_distance = 0.25;
    tree_ggd.name = "single_tree";

    Tree single_tree;
    gen.create_grove(tree_ggd, &single_tree, *data.heightmap);
    packer.pack_grove(ggd, res, *debugVisualizer, &single_tree, data.heightmap, visualize_voxels);

    print_alloc_info();
    distibutionGenerator.d = nullptr;
    dd.clear();
}
void generate_grove_renderer()
{
  std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
  if (pres == GroveRenderer::Precision::LOW)
    LODs_dists.back() = -10;
  data.groveRenderer = new GroveRenderer(&grove, &ggd, 5, LODs_dists, print_perf, pres);
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

  srand(time(NULL));

  data.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(2000, 2000), 5);
  data.heightmap->random_generate(0, 1, 50);

  if (generation_needed)
  {
    ggd = config.get_ggd(grove_type_name);
    
    if (parameter_selection)
    {
      std::function<void(TreeStructureParameters &, GrovePacked &)> _generate = generate_single_tree;
      ParameterSelector sel(_generate);
      TreeStructureParameters &start = ggd.types[0].params;
      Quality imp_qual = ggd.impostor_quality;
      ggd.impostor_quality = Quality::ULTRALOW;
      //ggd.bill_1_quality = Quality::ULTRALOW;
      //ggd.bill_2_quality = Quality::ULTRALOW;
      sel.select(start,SelectionType::SimulatedAnnealing, MetricType::ImpostorSimilarity);
      ggd.impostor_quality = imp_qual;
    }
    
    if (statistics_run)
    {
      ggd.clustering_max_individual_distance = statRunLaunchParams.max_ind_dist;
      ggd.trees_count = statRunLaunchParams.trees;
    }

    generate_grove();

    if (visualize_initial_voxels)
    {
      debugVisualizer->visualize_light_voxels(gen.voxels);
    }
  }
  if (saving_needed)
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
  generate_grove_renderer();

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
  textureManager.save_bmp(grove.impostors[1].atlas.tex(0),"impostor");

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

  if (false && need_initialization && appContext.renderMode != RenderMode::Rendering)
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
      generate_grove_renderer();
      appContext.regeneration_needed = false;
    }
    if (appContext.renderMode == RenderMode::StartingScreen)
      start_screen_pipeline();
    else if (appContext.renderMode == RenderMode::Rendering)
      simple_render_pipeline();
  };

  Tiny::view.interface = [&]()
  {
    bool p = true;
    ImGui::ShowDemoWindow(&p);
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
