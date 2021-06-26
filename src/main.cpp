
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

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
AppContext appContext;
ShadowMap shadowMap;
DefferedTarget defferedTarget;

Tree t[MAX_TREES];
mygen::Tree ttt;
mygen::TreeGenerator gen(ttt);
DebugVisualizer *debugVisualizer = nullptr;
GrovePacked grove;
GroveGenerationData ggd;
GroveRenderer *GR = nullptr;

int main(int argc, char *argv[])
{
  bool generation_needed = false;
  bool saving_needed = false;
  bool loading_needed = false;
  bool print_perf = false;
  bool only_gen = false;
  bool visualize_voxels = false;
  bool statistics_run = false;
  struct StatRunLaunchParams
  {
    int trees = 1;
    float max_ind_dist = 0.7;
  } statRunLaunchParams;
  GroveRenderer::Precision pres = GroveRenderer::MEDIUM;
  std::string grove_type_name = "default";
  std::string save_path = ".";
  std::string load_path = ".";
  int k = 1;
  while (k<argc)
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
        grove_type_name = argv[k+1];
      k+=2;
    }
    else if (std::string(argv[k]) == "-perf")
    {
      print_perf = true;
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_voxels")
    {
      visualize_voxels = true;
      k++;
    }
    else if (std::string(argv[k]) == "-only_gen")
    {
      only_gen = true;
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
        save_path = argv[k+1];
      k+=2;
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
        load_path = argv[k+1];
      k+=2;
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
        statRunLaunchParams.trees = atoi(argv[k+1]);
        statRunLaunchParams.max_ind_dist = atof(argv[k+2]);
        statistics_run = true;
        //only_gen = true;
        generation_needed = true;
        debug("stat run %d trees %f MID\n",statRunLaunchParams.trees,statRunLaunchParams.max_ind_dist);
        k+=3;
      }
    }
    else
    {
      logerr("unknown command \"%s\". Write -h for help",argv[k]);
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
    logerr("You choosed to do nothing. Program will be closed");
    return 0;
  }
  //base initialization
  glewInit();
  Tiny::view.lineWidth = 1.0f;
  Tiny::window("Procedural Tree", appContext.WIDTH, appContext.HEIGHT);
  Tiny::event.handler = [&](){eventHandler(appContext,Tiny::event);};
  textureManager = TextureManager("./resources/textures/");

  config.load_config();
  config.load_ggds();

  appContext.camera.pos = glm::vec3(-300,70,0);
  appContext.light.dir = glm::normalize(glm::vec3(-0.2,0.5,-0));
  appContext.light.color = glm::vec3(0.99,0.9,0.7);
  appContext.light.intensity = 1;
  appContext.light.ambient_q = 0.37;
  appContext.light.diffuse_q = 0.63;
  appContext.light.specular_q = 0.0;
  appContext.light.has_shadow_map = true;
  appContext.light.shadow_map_size = glm::vec2(4096, 4096);
  shadowMap.create(4096,4096);
  defferedTarget.create(2*1920,2*1080);
  defferedTarget.set_clear_color(glm::vec4(0.0,0.0,0.0,0.0));
  HBAORenderer hbaoRenderer;
  hbaoRenderer.create(800,600);
  Cubemap cubemap = Cubemap(1920, 1080);

  PostFx defferedLight = PostFx("deffered_light.fs");
  Shader defaultShader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  Shader depth({"depth.vs", "depth.fs"}, {"in_Position"});
  Shader debugShader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), &defaultShader);
  srand(time(NULL));
  std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
  if (pres == GroveRenderer::Precision::LOW)
    LODs_dists.back() = -10;
  Heightmap h = Heightmap(glm::vec3(0,0,0),glm::vec2(2000,2000),5);
  h.random_generate(0,1,50);
  GrovePacker packer;


  if (generation_needed)
  {
    ggd = config.get_ggd(grove_type_name);
    
    if (statistics_run)
    {
      ggd.clustering_max_individual_distance = statRunLaunchParams.max_ind_dist;
      ggd.trees_count = statRunLaunchParams.trees;
    }

    gen.create_grove(ggd, t, h);
    //Proctree::create_grove(ggd,t,h);
    packer.pack_grove(ggd,grove,*debugVisualizer, t,&h, visualize_voxels);
    distibutionGenerator.clear();
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
          logerr("path %s represents existing file. Can not save here",save_path.c_str());
          status = false;
        }
      } 
      if (status)
      {
        boost::filesystem::create_directory(save_path);
        boost::filesystem::permissions(save_path,boost::filesystem::perms::all_all);
      }
    }
    catch(const std::exception& e)
    {
      status = false;
      std::cerr << e.what() << '\n';
    }
    if (status)
    {
      std::string f_path = save_path + std::string("/grove.dat");
      FILE *f = fopen(f_path.c_str(), "wb");
      saver::set_textures_path(save_path);
      saver::save(f,grove);
      fclose(f);
    }
    else
    {
      logerr("error occured while saving to path %s",save_path.c_str());
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
          logerr("given load path %s is not a directory",load_path.c_str());
        }
      }
      else
      {
        logerr("given load path %s does not exist",load_path.c_str());
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
    }
  }
  GroveRenderer groveRenderer = GroveRenderer(&grove, &ggd, 5, LODs_dists, print_perf, pres);
  GR = &groveRenderer;
  for (int i=0;i<ggd.obstacles.size();i++)
    debugVisualizer->add_bodies(ggd.obstacles[i],1);
  TerrainRenderer tr = TerrainRenderer(h,glm::vec3(0,0,0),glm::vec2(2500,2500),glm::vec2(25,25));
  
  gltf::GeneralGltfWriter ggw;

  ggw.add_model(debugVisualizer->debugModels[0]);
  ggw.add_model(debugVisualizer->debugModels[1]);
  ggw.add_model(tr.flat_terrain);
  ggw.convert_to_gltf("terrain");
  ggw.clear();

  HeightmapTex ht = HeightmapTex(h,2048,2048);
  GrassRenderer gr = GrassRenderer();
  std::chrono::steady_clock::time_point t1, t_prev = std::chrono::steady_clock::now();
  float mu = 0.99;
  int frame = 0;
  bool regenerate_shadows = true;
  Timestamp ts;
  if (only_gen)
  {
    debug("Grove successfully generated. Exiting.");
    return 0;
  }
  Tiny::view.pipeline = [&]() 
  {
    auto &ctx = appContext;

    ctx.fpsCounter.tick();
    if (ctx.fpsCounter.get_frame_n() % 100 == 0 && print_perf)
    {
      fprintf(stderr,"FPS: %4.1f\n",ctx.fpsCounter.get_average_fps());
    }
    if (regenerate_shadows)
    {
      //shadows pass
      regenerate_shadows = false;
      shadowMap.use(ctx.light);
      glm::mat4 sh_viewproj = shadowMap.get_transform();

      groveRenderer.render(groveRenderer.get_max_LOD(), shadowMap.get_projection(), 
                          shadowMap.get_view(),ctx.camera,
                          glm::vec2(shadowMap.SHADOW_WIDTH,shadowMap.SHADOW_HEIGHT), 
                          ctx.light, ctx.groveRendererDebugParams,sh_viewproj,0,true);
      tr.render(shadowMap.get_projection(), shadowMap.get_view(),shadowMap.get_transform(),0,
                ctx.camera.pos,ctx.light, true);
      gr.render(shadowMap.get_projection(), shadowMap.get_view(),shadowMap.get_transform(),0,
                ctx.camera.pos, ht, ctx.light, true);

      shadowMap.start_trans_pass();
      groveRenderer.render(groveRenderer.get_max_LOD(), shadowMap.get_projection(), 
                          shadowMap.get_view(),ctx.camera,
                          glm::vec2(shadowMap.SHADOW_WIDTH,shadowMap.SHADOW_HEIGHT), 
                          ctx.light, ctx.groveRendererDebugParams,sh_viewproj,0,true);
      gr.render(shadowMap.get_projection(), shadowMap.get_view(),shadowMap.get_transform(),0,
                ctx.camera.pos, ht, ctx.light, true);
      shadowMap.finish_trans_pass();

      shadowMap.blur();
    }
    //Tiny::view.target(glm::vec3(0.6, 0.7, 1));
    defferedTarget.target();

    if (ctx.render_mode <= ctx.DEBUG_RENDER_MODE)
    {
      debugShader.use();
      debugShader.uniform("projectionCamera", ctx.projection * ctx.camera.camera());
      if (ctx.render_mode == ctx.DEBUG_RENDER_MODE)
      {
        debugShader.texture("tex", textureManager.get(ctx.debug_tex));
        debugShader.texture("tex_arr", textureManager.get(ctx.debug_tex));
        debugShader.uniform("need_tex",true);
        debugShader.uniform("need_arr_tex",false);
        debugShader.uniform("need_coord",false);
        debugShader.uniform("slice",0);
      }
      else if (ctx.render_mode == ctx.ARRAY_TEX_DEBUG_RENDER_MODE)
      {
        debugShader.texture("tex", textureManager.get_arr(ctx.debug_tex));
        debugShader.texture("tex_arr", textureManager.get_arr(ctx.debug_tex));
        debugShader.uniform("need_tex",false);
        debugShader.uniform("need_arr_tex",true);
        debugShader.uniform("need_coord",false);
        debugShader.uniform("slice", ctx.debug_layer);
      }
    }
    else
    {
      //depth prepass
      /*tr.render(ctx.projection * ctx.camera.camera(),shadowMap.get_transform(),shadowMap.getTex(),
                ctx.camera.pos,ctx.light, true);

      glClearColor(clearcolor.x, clearcolor.y, clearcolor.z, 1.0f);*/
      //color pass

      tr.render(ctx.projection, ctx.camera.camera(),shadowMap.get_transform(),0*shadowMap.getTex(),
                ctx.camera.pos,ctx.light);
      gr.render(ctx.projection, ctx.camera.camera(),shadowMap.get_transform(),0*shadowMap.getTex(),
               ctx.camera.pos, ht, ctx.light);
      if (ctx.render_mode != 2)
      {
        groveRenderer.render(ctx.forced_LOD, ctx.projection, ctx.camera.camera(),ctx.camera,
                             glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT), ctx.light, 
                             ctx.groveRendererDebugParams, shadowMap.get_transform(),0*shadowMap.getTex());  
      }
      debugVisualizer->render(ctx.projection * ctx.camera.camera(),ctx.render_mode);
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

  cubemap.render(ctx.projection, ctx.camera.camera(),ctx.camera);
  glm::vec3 ads = glm::vec3(ctx.light.ambient_q,ctx.light.diffuse_q,ctx.light.specular_q);
  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  defferedLight.use();
    defferedLight.get_shader().texture("colorTex",defferedTarget.get_color());
    defferedLight.get_shader().texture("normalsTex",defferedTarget.get_normals());
    defferedLight.get_shader().texture("viewPosTex",defferedTarget.get_view_pos());
    defferedLight.get_shader().texture("worldPosTex",defferedTarget.get_world_pos());
    defferedLight.get_shader().texture("aoTex",hbaoRenderer.get_tex());
    defferedLight.get_shader().texture("shadowMap",shadowMap.getTex());
    defferedLight.get_shader().texture("cubeTex", cubemap.get_tex());
    defferedLight.get_shader().uniform("dir_to_sun",ctx.light.dir);
    defferedLight.get_shader().uniform("camera_pos",ctx.camera.pos);
    defferedLight.get_shader().uniform("ambient_diffuse_specular",ads);
    defferedLight.get_shader().uniform("light_color",ctx.light.color);
    defferedLight.get_shader().uniform("need_shadow",shadowMap.getTex() != 0);
    defferedLight.get_shader().uniform("shadow_mat",shadowMap.get_transform());
    defferedLight.get_shader().uniform("sts_inv",1.0f/ctx.light.shadow_map_size);
    defferedLight.render();
  };

  Tiny::loop([&]() {});
  {
    for (int i=0;i<ggd.obstacles.size();i++)
    {
      delete ggd.obstacles[i];
      ggd.obstacles[i] = nullptr;
    }
    ggd.obstacles.clear();
  }
  Tiny::quit();

  return 0;
}
