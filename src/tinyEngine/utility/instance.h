#pragma once
#include <vector>
#include "model.h"
#include <GL/glew.h>

class Instance{
public:

  Instance(Primitive* _m){
    m = _m;
  };

  ~Instance(){
    if (buffer)
      delete(buffer);
    for(auto& b: instances){
        glDeleteBuffers(1, &b);
    }
  }

  Primitive* m;                     //Instanced Render Model (must be derived from primitive)
  std::vector<GLuint> instances;    //Instance VBO Pointers
  unsigned int SIZE;                //Number of Instances

  std::vector<glm::mat4> *buffer = nullptr;
  void addBuffer(std::vector<glm::mat4>& buf);
  void addBufferCopy(std::vector<glm::mat4>& buf);
  void updateBuffer(std::vector<glm::mat4>& buf, int index);
  void configBuffer(GLuint instance, std::vector<glm::mat4>& buf);
  void render(GLenum mode = GL_TRIANGLE_STRIP); //Default because of primitive models
};