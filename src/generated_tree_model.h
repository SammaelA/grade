#pragma once
#include "tinyEngine/TinyEngine.h"
#include "generated_tree.h"
#include "glm/trigonometric.hpp"
#include "tinyEngine/helpers/image.h"
#include "tinyEngine/helpers/color.h"
#include "tinyEngine/helpers/helper.h"
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#define PI 3.14159265f

const int WIDTH = 1200;
const int HEIGHT = 800;

float zoom = 0.5;
float zoomInc = 0.99;
float rotation = 0.0f;
glm::vec2 mousePos = glm::vec2(-1, -1);
float yaw = 0, pitch = 0;
glm::vec3 cameraPos = glm::vec3(0, 10, 10);
glm::vec3 cameraFront = glm::vec3(0, 0, -10);
glm::vec3 up = glm::vec3(0, 1, 0);
glm::mat4 camera = glm::lookAt(cameraPos, cameraPos + cameraFront, up);
glm::mat4 projection = glm::perspective(glm::radians(90.0f), (float)WIDTH / HEIGHT, 1.0f, 3000.0f);

bool paused = false;
bool autorotate = true;
bool drawwire = false;
bool drawtree = true;
bool drawleaf = true;

float leafcolor[3] = {0.82, 0.13, 0.23};
float treecolor[3] = {60. / 255, 30. / 255, 0.0};
float wirecolor[3] = {0.00, 0.00, 0.00};
float backcolor[3] = {0.80, 0.80, 0.80};
float lightcolor[3] = {1.00, 1.00, 1.00};
float leafopacity = 0.9;
int leafmindepth = 8;
float treeopacity = 1.0;
float treescale[2] = {1.0f, 1.0f};

int leafcount = 10;
float leafsize = 5.0;
float taper = 0.6;
float leafspread[3] = {50.0, 50.0, 50.0};

float growthrate = 5.0;
float passratio = 0.3;
float splitdecay = 1E-2;
float directedness = 0.5;
int localdepth = 4;
bool conservearea = true;

bool drawshadow = true;
bool selfshadow = true;
bool leafshadow = true;
glm::vec3 lightpos = glm::vec3(200);
glm::mat4 bias = glm::mat4(
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.5, 0.0,
    0.5, 0.5, 0.5, 1.0);
glm::mat4 lproj = glm::ortho(-1000.0f, 1000.0f, -1000.0f, 1000.0f, -200.0f, 2000.0f);
glm::mat4 lview = glm::lookAt(lightpos, glm::vec3(0), glm::vec3(0, 1, 0));
Tree t[100];
TreeGenerator gen(t[0]);
BillboardCloudRaw *cloud = nullptr;

void setup()
{
  //projection = glm::ortho(-(float)Tiny::view.WIDTH*zoom, (float)Tiny::view.WIDTH*zoom, -(float)Tiny::view.HEIGHT*zoom, (float)Tiny::view.HEIGHT*zoom, -500.0f, 800.0f);
  srand(time(NULL));
  float bp[] = {0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10};
  for (int i = 0; i < 2; i++)
  {
    TreeStructureParameters par;
    t[i].pos += glm::vec3(100.0 * (i % 10), 0, 100.0 * (i / 10));
    /*int row = i % 10;
    int column = i / 10;
    switch(column)
    {
      case 0:
        par.base_branch_feed = 20*(i+1);
        par.base_seg_feed = 5*(i+1);
        
        break;
      case 1:
        par.feed_distribution_min_weight = 0.05 + 0.1*row;
        break;
      case 2:
        par.top_growth_bonus = 0.04*(i+1);
        break;
      case 3:
        par.base_branch_feed = 5000 + 1500*i;
        break;
      case 4:
        par.base_seg_feed = 300 + 100*i;
        break;

      case 5:
        par.seg_spread = 0.02 + 0.02*row;
        break;
      case 6:
        par.seg_phototrop = 0.02 + 0.02*row;
        break;
      case 7:
        par.seg_gravitrop = 0.02 + 0.02*row;
        break;
      case 8:
        par.seg_dir_conserv = 0.2 + 0.2*row;
        break;
      
      default:
        break;
    }*/
    t[i].id = i;
    gen.plant_tree(t[i], par);
  }
}

// Event Handler
std::function<void()> eventHandler = [&]() {
  /*if(Tiny::event.scroll.posy){
    zoom /= zoomInc;
    //projection = glm::ortho(-(float)Tiny::view.WIDTH*zoom, (float)Tiny::view.WIDTH*zoom, -(float)Tiny::view.HEIGHT*zoom, (float)Tiny::view.HEIGHT*zoom, -500.0f, 800.0f);
  }
  if(Tiny::event.scroll.negy){
    zoom *= zoomInc;
    //projection = glm::ortho(-(float)Tiny::view.WIDTH*zoom, (float)Tiny::view.WIDTH*zoom, -(float)Tiny::view.HEIGHT*zoom, (float)Tiny::view.HEIGHT*zoom, -500.0f, 800.0f);
  }
  if(Tiny::event.scroll.posx){
    rotation += 1.5f;
    if(rotation < 0.0) rotation = 360.0 + rotation;
    else if(rotation > 360.0) rotation = rotation - 360.0;
    //camera = glm::rotate(camera, glm::radians(1.5f), glm::vec3(0.0f, 1.0f, 0.0f));
  }
  if(Tiny::event.scroll.negx){
    rotation -= 1.5f;
    if(rotation < 0.0) rotation = 360.0 + rotation;
    else if(rotation > 360.0) rotation = rotation - 360.0;
    //camera = glm::rotate(camera, glm::radians(-1.5f), glm::vec3(0.0f, 1.0f, 0.0f));
  }*/
  float nx = Tiny::event.mouse.x;
  float ny = Tiny::event.mouse.y;
  GLfloat xoffset = nx - mousePos.x;
  GLfloat yoffset = mousePos.y - ny;

  GLfloat sensitivity = 0.1;
  xoffset *= sensitivity;
  yoffset *= sensitivity;

  yaw += xoffset;
  pitch += yoffset;

  if (pitch > 89.0f)
    pitch = 89.0f;
  if (pitch < -89.0f)
    pitch = -89.0f;
  glm::vec3 front;
  front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
  front.y = sin(glm::radians(pitch));
  front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
  cameraFront = glm::normalize(front);
  mousePos = glm::vec2(nx, ny);
  //Pause Toggle
  float speed = 2.0f;
  glm::vec3 cameraPerp = glm::normalize(glm::cross(cameraFront, up));
  if (Tiny::event.active[SDLK_w])
    cameraPos += speed * cameraFront;
  if (Tiny::event.active[SDLK_s])
    cameraPos -= speed * cameraFront;

  if (Tiny::event.active[SDLK_a])
    cameraPos -= speed * cameraPerp;
  if (Tiny::event.active[SDLK_d])
    cameraPos += speed * cameraPerp;

  if (Tiny::event.active[SDLK_e])
    cameraPos += speed * up;
  if (Tiny::event.active[SDLK_c])
    cameraPos -= speed * up;

  if (!Tiny::event.press.empty())
  {

    if (Tiny::event.press.back() == SDLK_p)
      paused = !paused;
    else if (Tiny::event.press.back() == SDLK_r)
      autorotate = !autorotate;

    //Regrow
    else if (Tiny::event.press.back() == SDLK_r)
    {
    }
  }
  camera = glm::lookAt(cameraPos, cameraPos + cameraFront, up);
};

std::function<void(Model *)> _ce = [&](Model *h) {};
std::function<void(Model *)> _construct = [&](Model *h) {
  for (int i = 1; i < 2; i++)
  {
    gen.grow_tree(t[i]);
    gen.tree_to_model(t[i], h, false, cloud);
  }
};
std::function<void(Model *)> _construct_leaves = [&](Model *h) {
  for (int i = 1; i < 2; i++)
  {
    gen.tree_to_model(t[i], h, true);
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
