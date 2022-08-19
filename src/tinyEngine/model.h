#pragma once
#include "tinyEngine/resources.h"
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <functional>
#include <iostream>
#include <chrono>
#include "common_utils/utility.h"
#include "texture.h"

struct Mesh
{
  std::vector<GLfloat>  positions;//vec3
  std::vector<GLfloat>  normals;//vec3
  std::vector<GLfloat>  colors;//vec4
  std::vector<GLfloat>  tangents;//vec4, optional
  std::vector<GLuint>   indices;//3*triangles_cnt
  std::vector<int>   mat_indicies;//triangles_cnt, optional
  bool indexed = true;
};

struct Material
{
  glm::vec3 Ka = glm::vec3(1, 1, 1); // Ambient Color
  glm::vec3 Kd = glm::vec3(1, 1, 1); // Diffuse Color
  glm::vec3 Ks = glm::vec3(1, 1, 1); // Specular Color
  float Ns = 32;                     // Specular Exponent
  float Ni = 1;                      // Optical Density
  float d = 1;                       // Dissolve
  int illum = 0;                     // Illumination Model

  Texture map_Ka;   // Ambient Texture Map
  Texture map_Kd;   // Diffuse Texture Map
  Texture map_Ks;   // Specular Texture Map
  Texture map_Ns;   // Specular Hightlight Map
  Texture map_d;    // Alpha Texture Map
  Texture map_bump; // Bump Map

  Material() {};
  Material(Texture &tex)
  {
    //simple material with only albedo texture
    map_Ka = tex;
  }
};

class Model: public Mesh //render-ready model that creates it's own vao and vbo's 
{
public:
  GLuint vao;
  std::vector<GLuint> vbo;
  GLuint ibo;
  size_t SIZE = 4;
  glm::mat4 model = glm::mat4(1.0);

  virtual int get_size() {return positions.size()/3;}

  void bindf(int index, int count, int size, float* data)
  {
    glBindBuffer(GL_ARRAY_BUFFER, vbo[index]);
    glBufferData(GL_ARRAY_BUFFER, count*sizeof(float), data, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(index);
    glVertexAttribPointer(index, size, GL_FLOAT, GL_FALSE, 0, 0);
  }
  Model()
  {
    vao = create_vertex_array();
    glBindVertexArray(vao);
    ibo = create_buffer();
    for(int i = 0; i < 3; i++)//positions, normals, colors
    {
      GLuint nvbo = create_buffer();
      vbo.push_back(nvbo);
    }
  }

  ~Model()
  {
    glDisableVertexAttribArray(vao);
    delete_vertex_array(vao);
    delete_buffer(ibo);
    for (auto vb : vbo)
      delete_buffer(vb);
  }

  void update()
  {
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    glBindVertexArray(vao);

    bindf(0, positions.size(), 3, &positions[0]);
    if(!normals.empty()) 
      bindf(1, normals.size(), 3, &normals[0]);
    if(!colors.empty()) 
      bindf(2, colors.size(), 4, &colors[0]);

    SIZE = indices.size();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, SIZE*sizeof(GLuint), &indices[0], GL_DYNAMIC_DRAW);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  }

  void construct(std::function<void(Model*)> constructor)
  {
    positions.clear();
    normals.clear();
    colors.clear();
    indices.clear();

    (constructor)(this);  //Call user-defined constructor
    update();                //Update VAO / VBO / IBO
  }

  void render(GLenum mode = GL_TRIANGLES)
  {
    glBindVertexArray(vao);
    if(indexed)
    {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
      glDrawElements(mode, SIZE, GL_UNSIGNED_INT, 0);
    }
    else 
      glDrawArrays(mode, 0, positions.size()/3);
  }
};

struct ComplexModel //an object that contains different models with different material
{
  std::vector<Model *> models;
  std::vector<Material> materials;
  void update()
  {
    for (Model *m : models)
      m->update();
  }
};

//Primitive Shapes (Pre-Made)

struct Square2D: Model
{
  GLfloat vert[8] = {-1.0,  1.0, -1.0, -1.0,  1.0,  1.0,  1.0, -1.0};
  GLfloat tex [8] = { 0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  1.0,  1.0};

  Square2D():Model()
  {
    indexed = false;
    bindf(0, 8, 2, &vert[0]);
    bindf(1, 8, 2, &tex[0]);
  }
};

struct Square3D: Model
{
  GLfloat vert[12] =  {-1.0,  1.0,  0.0, -1.0, -1.0,  0.0,  1.0,  1.0,  0.0,  1.0, -1.0,  0.0};
  GLfloat tex [8]   = { 0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  1.0,  1.0};

  Square3D():Model()
  {
    indexed = false;
    bindf(0, 12, 3, &vert[0]);
    bindf(1, 8,  2, &tex[0]);
  }
};

struct Cube: Model
{
  GLfloat vert[3*6*6] = {
    // positions          
    -1.0f,  1.0f, -1.0f,
    -1.0f, -1.0f, -1.0f,
     1.0f, -1.0f, -1.0f,
     1.0f, -1.0f, -1.0f,
     1.0f,  1.0f, -1.0f,
    -1.0f,  1.0f, -1.0f,

    -1.0f, -1.0f,  1.0f,
    -1.0f, -1.0f, -1.0f,
    -1.0f,  1.0f, -1.0f,
    -1.0f,  1.0f, -1.0f,
    -1.0f,  1.0f,  1.0f,
    -1.0f, -1.0f,  1.0f,

     1.0f, -1.0f, -1.0f,
     1.0f, -1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
     1.0f,  1.0f, -1.0f,
     1.0f, -1.0f, -1.0f,

    -1.0f, -1.0f,  1.0f,
    -1.0f,  1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
     1.0f, -1.0f,  1.0f,
    -1.0f, -1.0f,  1.0f,

    -1.0f,  1.0f, -1.0f,
     1.0f,  1.0f, -1.0f,
     1.0f,  1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
    -1.0f,  1.0f,  1.0f,
    -1.0f,  1.0f, -1.0f,

    -1.0f, -1.0f, -1.0f,
    -1.0f, -1.0f,  1.0f,
     1.0f, -1.0f, -1.0f,
     1.0f, -1.0f, -1.0f,
    -1.0f, -1.0f,  1.0f,
     1.0f, -1.0f,  1.0f
};
  GLfloat tex [8]   = { 0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  1.0,  1.0};

  Cube():Model()
  {
    indexed = false;
    bindf(0, 3*6*6, 3, &vert[0]);
    bindf(1, 8,  2, &tex[0]);
    SIZE = 36;
  }
};