#pragma once
#include <vector>
#include "model.h"
#include <GL/glew.h>

class Instance : Countable{
public:

  Instance(Primitive* _m);
  ~Instance();

  Primitive* m;                     //Instanced Render Model (must be derived from primitive)
  std::vector<GLuint> instances;    //Instance VBO Pointers
  unsigned int SIZE;                //Number of Instances

  std::vector<std::vector<glm::mat4> *> buffers;
  std::vector<std::vector<glm::vec4> *> vectors;

  void addBuffer(std::vector<glm::mat4>& buf);
  void addBufferCopy(std::vector<glm::mat4>& buf);
  void updateBuffer(std::vector<glm::mat4>& buf, int index);
  void configBuffer(GLuint instance, std::vector<glm::mat4>& buf);

  void addBuffer(std::vector<glm::vec4>& buf);
  void addBufferCopy(std::vector<glm::vec4>& buf);
  void updateBuffer(std::vector<glm::vec4>& buf, int index);
  void configBuffer(GLuint instance, std::vector<glm::vec4>& buf);

  void render(GLenum mode = GL_TRIANGLE_STRIP); //Default because of primitive models
};