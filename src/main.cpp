
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

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Camera camera;

const int WIDTH = 1200;
const int HEIGHT = 800;
int treecount = 0;
int cloudnum = 1;
bool draw_clusterized = true;
int cur_tree = 0;
Tree t[101];
TreeGenerator gen(t[100]);
DebugVisualizer debugVisualizer;
GrovePacked grove;
BillboardCloudRaw::RenderMode mode = BillboardCloudRaw::ONLY_SINGLE;
glm::vec2 mousePos = glm::vec2(-1, -1);
glm::mat4 projection = glm::perspective(glm::radians(90.0f), (float)WIDTH / HEIGHT, 1.0f, 3000.0f);

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
  treecount = 2;
  TreeStructureParameters params;
  srand(time(NULL));
  float bp[] = {0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10};
  gen.create_grove(params, treecount, grove, debugVisualizer, t);
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

  if (Tiny::event.active[SDLK_u])
    cloudnum = 0;
  if (Tiny::event.active[SDLK_i])
    cloudnum = 1;
  if (Tiny::event.active[SDLK_o])
    cloudnum = 2;
  if (Tiny::event.active[SDLK_p])
    cloudnum = 3;
  if (Tiny::event.active[SDLK_k])
  {
    draw_clusterized = true;
  }
  if (Tiny::event.active[SDLK_l])
  {
    draw_clusterized = false;
  }
  if (!Tiny::event.press.empty())
  {
  }
};

std::function<void(Model *m)> construct_floor = [&](Model *h) {
  float floor[12] = {
      -100.0,
      10.0,
      -100.0,
      -100.0,
      0.0,
      100.0,
      100.0,
      10.0,
      -100.0,
      100.0,
      0.0,
      100.0,
  };

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

    h->colors.push_back(floorcolor.x);
    h->colors.push_back(floorcolor.y);
    h->colors.push_back(floorcolor.z);
    h->colors.push_back(1.0);
  }
};

int main(int argc, char *args[])
{
  glewInit();
  Tiny::view.lineWidth = 1.0f;
  Tiny::window("Procedural Tree", WIDTH, HEIGHT);
  Tiny::event.handler = eventHandler;
  textureManager = TextureManager("resources/textures/");
  Texture &tex = textureManager.get("woodd");
  Texture &wood = textureManager.get("wood");
  TreeStructureParameters par;
  for (int i = 0; i < 100; i++)
  {
    t[i] = Tree();
    t[i].leaf = &tex;
    t[i].wood = &wood;
    t[i].params = par;
  }

  Model floor(construct_floor);

  Shader particleShader({"particle.vs", "particle.fs"}, {"in_Quad", "in_Tex", "in_Model"});
  Shader defaultShader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  Shader depth({"depth.vs", "depth.fs"}, {"in_Position"});
  BillboardTiny shadow(1600, 1600, false);
  debugVisualizer = DebugVisualizer(&wood, &defaultShader);
  setup();
  GroveRenderer groveRenderer = GroveRenderer(&grove, 2);
  groveRenderer.pwood = &wood;
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
    defaultShader.use();
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
      groveRenderer.render(cloudnum, projection * camera.camera());
    }
    else
    {
      defaultShader.uniform("drawfloor", false);
      defaultShader.texture("tex", wood);
      defaultShader.uniform("wireframe", false);
      glm::mat4 prc = projection * camera.camera();
      for (int i = 0; i < treecount; i++)
      {
        t[i].render(defaultShader, cloudnum, prc);
      }
      debugVisualizer.render(prc);
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
  if (wood)
    defaultShader.texture("tex", *wood);
  defaultShader.uniform("model", models[cloudnum]->model);
  models[cloudnum]->update();
  models[cloudnum]->render(GL_TRIANGLES);
  if (cloudnum == 3 && models.size() >= 5)
  {
    if (leaf)
      defaultShader.texture("tex", *leaf);
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