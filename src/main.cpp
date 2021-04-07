
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
#include "camera.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "tinyEngine/utility.h"
#include "grove.h"
#include "tinyEngine/save_utils/config.h"
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include "terrain.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
Camera camera;

const int WIDTH = 1200;
const int HEIGHT = 800;
float avg_fps = 100;
int treecount = 0;
int cloudnum = -1;
bool draw_clusterized = true;
int cur_tree = 0;
const int DEBUG_RENDER_MODE = -2;
const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
int render_mode = -1;
int debug_tex = 0;
int debug_layer = 0;
int debug_model_focus = 0;
bool debug_need_focus = true;
float rot_a;
Tree t[MAX_TREES];
TreeGenerator gen(t[100]);
DebugVisualizer *debugVisualizer = nullptr;
GrovePacked grove;
GroveGenerationData ggd;
BillboardCloudRaw::RenderMode mode = BillboardCloudRaw::ONLY_SINGLE;
glm::vec2 mousePos = glm::vec2(-1, -1);
glm::mat4 projection = glm::perspective(glm::radians(90.0f), (float)WIDTH / HEIGHT, 1.0f, 3000.0f);
GroveRenderer *GR = nullptr;
bool drawshadow = true;
glm::vec3 lightpos = glm::vec3(200);
glm::mat4 bias = glm::mat4(
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.5, 0.0,
    0.5, 0.5, 0.5, 1.0);
glm::mat4 lproj = glm::ortho(-1000.0f, 1000.0f, -1000.0f, 1000.0f, -200.0f, 2000.0f);
glm::mat4 lview = glm::lookAt(lightpos, glm::vec3(0), glm::vec3(0, 1, 0));


// Event Handler
std::function<void()> eventHandler = []() {
  float nx = Tiny::event.mouse.x;
  float ny = Tiny::event.mouse.y;
  GLfloat xoffset = nx - mousePos.x;
  GLfloat yoffset = mousePos.y - ny;

  GLfloat sensitivity = 0.1;
  xoffset *= sensitivity;
  yoffset *= sensitivity;

  camera.yaw += xoffset;
  camera.pitch += yoffset;

  if (camera.pitch > 89.0f)
    camera.pitch = 89.0f;
  if (camera.pitch < -89.0f)
    camera.pitch = -89.0f;
  glm::vec3 front;
  front.x = cos(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
  front.y = sin(glm::radians(camera.pitch));
  front.z = sin(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
  camera.front = glm::normalize(front);
  mousePos = glm::vec2(nx, ny);
  //Pause Toggle
  float speed = 0.2;
  glm::vec3 cameraPerp = glm::normalize(glm::cross(camera.front, camera.up));
  if (Tiny::event.active[SDLK_l])
  {
    camera = Camera();
    camera.pos = glm::vec3(-200,70,0);
  }
  if (Tiny::event.active[SDLK_w])
    camera.pos += speed * camera.front;
  if (Tiny::event.active[SDLK_s])
    camera.pos -= speed * camera.front;

  if (Tiny::event.active[SDLK_a])
    camera.pos -= speed * cameraPerp;
  if (Tiny::event.active[SDLK_d])
    camera.pos += speed * cameraPerp;

  if (Tiny::event.active[SDLK_e])
    camera.pos += speed * camera.up;
  if (Tiny::event.active[SDLK_c])
    camera.pos -= speed * camera.up;
  if (Tiny::event.active[SDLK_y])
    cloudnum = -1;
  if (Tiny::event.active[SDLK_u])
    cloudnum = 0;
  if (Tiny::event.active[SDLK_i])
    cloudnum = 1;
  if (Tiny::event.active[SDLK_o])
    cloudnum = 2;
  if (Tiny::event.active[SDLK_p])
    cloudnum = 3;
  if (Tiny::event.active[SDLK_LEFTBRACKET])
    cloudnum = 4;
  if (Tiny::event.active[SDLK_RIGHTBRACKET])
    cloudnum = 5;
  if (Tiny::event.active[SDLK_k])
  {
    draw_clusterized = true;
  }
  if (Tiny::event.active[SDLK_l])
  {
    //draw_clusterized = false;
  }
  if (Tiny::event.active[SDLK_r])
  {
    for (int i=0;i<MAX_TREES;i++)
    {
      t[i] = Tree();
    }
    if (GR)
    {
      GR->~GroveRenderer();
      GR = nullptr;
    }
  }
  if (Tiny::event.active[SDLK_m])
  {
    render_mode++;
    if (render_mode > DebugVisualizer::MAX_RENDER_MODE)
      render_mode = ARRAY_TEX_DEBUG_RENDER_MODE;
    Tiny::event.active[SDLK_m] = false;
    logerr("render mode %d",render_mode);
  }
  if (Tiny::event.active[SDLK_n])
  {
    if (render_mode <= DEBUG_RENDER_MODE)
      debug_tex++;
    Tiny::event.active[SDLK_n] = false;
  }
  if (Tiny::event.active[SDLK_b])
  {
    if (debug_need_focus)
      debug_model_focus++;
    if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
      debug_layer++;
    Tiny::event.active[SDLK_b] = false;
  }
    if (Tiny::event.active[SDLK_v])
  {
    if (debug_need_focus)
      debug_model_focus--;
    if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
      debug_layer--;
    Tiny::event.active[SDLK_v] = false;
  }
  if (Tiny::event.active[SDLK_f])
  {
    debug_need_focus = !debug_need_focus;
    render_mode = -1;
    Tiny::event.active[SDLK_f] = false;
  }
  if (!Tiny::event.press.empty())
  {

  }
};

std::function<void(Model *m)> construct_floor = [](Model *h) {
  float floor[12] = {
      -100.0,
      0.0,
      -100.0,
      -100.0,
      0.0,
      100.0,
      100.0,
      0.0,
      -100.0,
      100.0,
      0.0,
      100.0,
  };
  float colors[8] = {0,0,0,1,1,0,1,1};
  for (int i = 0; i < 12; i++)
    h->positions.push_back(floor[i]);

  h->indices.push_back(0);
  h->indices.push_back(1);
  h->indices.push_back(2);

  h->indices.push_back(1);
  h->indices.push_back(3);
  h->indices.push_back(2);

  glm::vec3 floorcolor = glm::vec3(0.1, 0.3, 0.1);

  for (int i = 0; i < 4; i++)
  {
    h->normals.push_back(0.0);
    h->normals.push_back(1.0);
    h->normals.push_back(0.0);

    h->colors.push_back(colors[2*i]);
    h->colors.push_back(colors[2*i + 1]);
    h->colors.push_back(0.0);
    h->colors.push_back(1.0);
  }
};

int main(int argc, char *argv[])
{
  bool generation_needed = false;
  bool saving_needed = false;
  bool loading_needed = false;
  bool print_perf = false;
  bool only_gen = false;
  bool visualize_voxels = false;
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
  Tiny::window("Procedural Tree", WIDTH, HEIGHT);
  Tiny::event.handler = eventHandler;
  textureManager = TextureManager("/home/sammael/study/bit_bucket/grade/resources/textures/");

  config.load_config();
  config.load_ggds();

  camera.pos = glm::vec3(-300,70,0);

  Model floor(construct_floor);

  Shader particleShader({"particle.vs", "particle.fs"}, {"in_Quad", "in_Tex", "in_Model"});
  Shader defaultShader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  Shader depth({"depth.vs", "depth.fs"}, {"in_Position"});
  Shader debugShader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  BillboardTiny shadow(1600, 1600, false);
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), &defaultShader);
  srand(time(NULL));
  std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
  if (pres == GroveRenderer::Precision::LOW)
    LODs_dists.back() = -10;
  Heightmap h = Heightmap(glm::vec3(0,0,0),glm::vec2(1000,1000),5);
  h.random_generate(0,1,50);
  if (generation_needed)
  {
    ggd = config.get_ggd(grove_type_name);
    gen.create_grove(ggd, grove, *debugVisualizer, t, &h, visualize_voxels);
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
  
  std::chrono::steady_clock::time_point t1, t_prev = std::chrono::steady_clock::now();
  float mu = 0.99;
  int frame = 0;
  Timestamp ts;
  if (only_gen)
  {
    logerr("Grove successfully generated. Exiting.");
    return 0;
  }
  Tiny::view.pipeline = [&]() {
    t_prev = t1;
    t1 = std::chrono::steady_clock::now();
    float frame_time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t_prev).count();
    frame_time = MAX(frame_time,0.1);
    avg_fps = mu*avg_fps + (1 - mu)*(1000/frame_time);
    frame++;
    if (frame % 100 == 0 && print_perf)
    {
      fprintf(stderr,"FPS: %4.1f\n",avg_fps);
    }
    rot_a += 0.01;
    //camera.pos = glm::vec3(250*sin(rot_a),100,250*cos(rot_a));
    //camera.front = glm::normalize(-camera.pos);

    Tiny::view.target(glm::vec3(0.6, 0.7, 1));
    if (render_mode <= DEBUG_RENDER_MODE)
    {
      debugShader.use();
      debugShader.uniform("projectionCamera", projection * camera.camera());
      if (render_mode == DEBUG_RENDER_MODE)
      {
        debugShader.texture("tex", textureManager.get(debug_tex));
        debugShader.texture("tex_arr", textureManager.get(debug_tex));
        debugShader.uniform("need_tex",true);
        debugShader.uniform("need_arr_tex",false);
        debugShader.uniform("need_coord",false);
        debugShader.uniform("slice",0);
      }
      else if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
      {
        debugShader.texture("tex", textureManager.get_arr(debug_tex));
        debugShader.texture("tex_arr", textureManager.get_arr(debug_tex));
        debugShader.uniform("need_tex",false);
        debugShader.uniform("need_arr_tex",true);
        debugShader.uniform("need_coord",false);
        debugShader.uniform("slice", debug_layer);
      }
      debugShader.uniform("model", floor.model);
      floor.render(GL_TRIANGLES);
    }
    else
    {
      tr.render(projection * camera.camera());
      if (draw_clusterized && render_mode != 2)
      {
        GroveRendererDebugParams dbgpar;
        dbgpar.need_focus_model = debug_need_focus;
        dbgpar.model_focused = debug_model_focus;
        groveRenderer.render(cloudnum, projection * camera.camera(),camera, 
        glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT), dbgpar);
        
      }
      debugVisualizer->render(projection * camera.camera(),render_mode);
    }
  };

  //Loop over Stuff
  Tiny::loop([&]() { /* ... */
                     floor.construct(construct_floor);
  });

  Tiny::quit();

  return 0;
}
void Tree::render(Shader &defaultShader, int cloudnum, glm::mat4 prc)
{
  if (models.size() == 0 || billboardClouds.size() == 0)
  {
    logerr("wtf empty tree id =  %d  %d %d\n", id, cloudnum, models.size(), billboardClouds.size());
    return;
  }

  if (cloudnum < 0)
    cloudnum = 0;

  else if (cloudnum >= billboardClouds.size())
    cloudnum = billboardClouds.size() - 1;

  defaultShader.use();
  defaultShader.texture("tex", wood);
  defaultShader.uniform("model", models[cloudnum]->model);
  models[cloudnum]->update();
  models[cloudnum]->render(GL_TRIANGLES);
  if (cloudnum == 3 && models.size() >= 5)
  {
    defaultShader.texture("tex", leaf);
    defaultShader.uniform("model", models[cloudnum + 1]->model);
    models[cloudnum + 1]->update();
    models[cloudnum + 1]->render(GL_TRIANGLES);
  }
  else
  {
    billboardClouds[cloudnum]->set_render_mode(mode);
    billboardClouds[cloudnum]->render(prc);
  }
}