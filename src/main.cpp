
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

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
Camera camera;

const int WIDTH = 1200;
const int HEIGHT = 800;
int treecount = 0;
int cloudnum = -1;
bool draw_clusterized = true;
int cur_tree = 0;
const int DEBUG_RENDER_MODE = -2;
const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
int render_mode = -1;
int debug_tex = 0;
int debug_layer = 0;
Tree t[101];
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

void setup()
{
  treecount = 3;
  TreeStructureParameters params1,params2,params3;
  params2.seg_len_mult = Parameter<float>(4, std::vector<float>{0.1, 1.75, 1, 0.55, 0.4});
  params2.base_seg_feed = Parameter<float>(100, std::vector<float>{20, 20, 120, 40, 30}, REGENERATE_ON_GET, new Uniform(-0, 0));
  params2.base_branch_feed = Parameter<float>(300, std::vector<float>{20, 30, 200, 40, 40}, REGENERATE_ON_GET, new Uniform(-00, 00));
  params2.min_branching_chance = Parameter<float>(0, std::vector<float>{1, 0.4, 0.5, 0.75, 0.7});
  params1 = config.get("default");
  params2 = config.get("bush");
  TreeTypeData ttd1(0,params1,std::string("wood"),std::string("leaf"));
  TreeTypeData ttd2(1,params2,std::string("wood2"),std::string("leaf2"));
  TreeTypeData ttd3(2,params3,std::string("wood3"),std::string("leaf2"));
  float r = sqrt(treecount) + 1;
  ggd.size = glm::vec3(40.0f * r + 250, 120, 40.0f * r + 250);
  ggd.pos = glm::vec3(0,0,0);
  ggd.trees_count = treecount;
  ggd.types = {ttd1,ttd2,ttd3};
  srand(time(NULL));
  gen.create_grove(config.get_ggd("bush_grove"), grove, *debugVisualizer, t);
}

// Event Handler
std::function<void()> eventHandler = [&]() {
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
  float speed = 2.0f;
  glm::vec3 cameraPerp = glm::normalize(glm::cross(camera.front, camera.up));
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
  if (Tiny::event.active[SDLK_k])
  {
    draw_clusterized = true;
  }
  if (Tiny::event.active[SDLK_l])
  {
    draw_clusterized = false;
  }
  if (Tiny::event.active[SDLK_r])
  {
    for (int i=0;i<101;i++)
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
    if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
      debug_layer++;
    Tiny::event.active[SDLK_b] = false;
  }
    if (Tiny::event.active[SDLK_v])
  {
    if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
      debug_layer--;
    Tiny::event.active[SDLK_v] = false;
  }
  if (!Tiny::event.press.empty())
  {

  }
};

std::function<void(Model *m)> construct_floor = [&](Model *h) {
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

  camera.pos = glm::vec3(-200,50,0);

  Model floor(construct_floor);

  Shader particleShader({"particle.vs", "particle.fs"}, {"in_Quad", "in_Tex", "in_Model"});
  Shader defaultShader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  Shader depth({"depth.vs", "depth.fs"}, {"in_Position"});
  Shader debugShader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  BillboardTiny shadow(1600, 1600, false);
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), &defaultShader);
  srand(time(NULL));
  std::vector<float> LODs_dists = {1000,500,300,100};

  if (generation_needed)
  {
    ggd = config.get_ggd(grove_type_name);
    gen.create_grove(ggd, grove, *debugVisualizer, t);
  }
  if (saving_needed)
  {
    struct stat sb;
    if (stat(save_path.c_str(), &sb) != 0)
    {
      mkdir(save_path.c_str(),0777);
    }
    if (!S_ISDIR(sb.st_mode))
    {
      logerr("given save path \"%s\" is not a directory",save_path.c_str());
    }
    std::string f_path = save_path + std::string("/grove.dat");
    FILE *f = fopen(f_path.c_str(), "wb");
    saver::set_textures_path(save_path);
    saver::save(f,grove);
    fclose(f);
  }
  if (loading_needed)
  {
    grove = GrovePacked();
    struct stat sb;
    if (stat(load_path.c_str(), &sb) != 0)
    {
      mkdir(load_path.c_str(),0777);
    }
    if (!S_ISDIR(sb.st_mode))
    {
      logerr("given load path \"%s\" is not a directory",load_path.c_str());
    }
    std::string f_path = load_path + std::string("/grove.dat");
    FILE *f = fopen(f_path.c_str(), "rb");
    saver::set_textures_path(load_path);
    saver::load(f,grove);
    ggd = config.get_ggd(grove.ggd_name);
    fclose(f);
  }
  GroveRenderer groveRenderer = GroveRenderer(&grove, &ggd, 4, LODs_dists);
  GR = &groveRenderer;




  Tiny::view.pipeline = [&]() {
    shadow.target();
    if (drawshadow)
    {
      depth.use();
      depth.uniform("dvp", lproj * lview);
      //	for (int i=0;i<0*treecount;i++)
      //	{
      //		t[i].render(defaultShader,cloudnum,lproj*lview);
      //	}
    }
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
      defaultShader.use();
      defaultShader.uniform("mult", 0);
      defaultShader.uniform("projectionCamera", projection * camera.camera());
      defaultShader.uniform("lightcolor", glm::vec3(1, 1, 1));
      defaultShader.uniform("lookDir", camera.front);
      defaultShader.uniform("lightDir", lightpos);

      defaultShader.uniform("drawshadow", drawshadow);
      if (drawshadow)
      {
        defaultShader.uniform("dbvp", bias * lproj * lview);
        defaultShader.texture("shadowMap", shadow.depth);
        defaultShader.uniform("light", lightpos);
      }

      defaultShader.uniform("drawfloor", true);
      defaultShader.uniform("drawcolor", glm::vec4(floor.colors[0], floor.colors[1], floor.colors[2], 1));
      defaultShader.uniform("model", floor.model);
      floor.render(GL_TRIANGLES);
      if (draw_clusterized)
      {
        groveRenderer.render(cloudnum, projection * camera.camera(),camera.pos,glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT));
      }
      else
      {//deprecated
      /*
        defaultShader.uniform("drawfloor", false);
        defaultShader.texture("tex", wood);
        defaultShader.uniform("wireframe", false);
        glm::mat4 prc = projection * camera.camera();
        for (int i = 0; i < treecount; i++)
        {
          t[i].render(defaultShader, cloudnum, prc);
        }
        debugVisualizer->render(prc,render_mode);*/
      }
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