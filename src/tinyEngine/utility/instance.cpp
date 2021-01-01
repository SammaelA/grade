#include "instance.h"

void Instance::addBuffer(std::vector<glm::mat4>& buf){
  GLuint instance;
  glGenBuffers(1, &instance);
  SIZE = buf.size();              //Update the Number of Instances

  glBindVertexArray(m->vao);
  glBindBuffer(GL_ARRAY_BUFFER, instance);  //Bind Instance Buffer and Data
  glBufferData(GL_ARRAY_BUFFER, SIZE*sizeof(glm::mat4), &buf[0], GL_STATIC_DRAW);

  configBuffer(instance,buf);
}

void Instance::configBuffer(GLuint instance, std::vector<glm::mat4>& buf){
  for(int i = 0; i < 4; i++)//For Matrices - Special Procedure
  {
    glEnableVertexAttribArray(m->vbo.size()+instances.size());
    glVertexAttribPointer(m->vbo.size()+instances.size(), 4, GL_FLOAT, GL_FALSE, 4*sizeof(glm::vec4), (void*)(i*sizeof(glm::vec4)));
    glVertexAttribDivisor(m->vbo.size()+instances.size(), 1);
    instances.push_back(instance);
  }
  /*for other types
  glEnableVertexAttribArray(m->vbo.size()+instances.size());
  glVertexAttribPointer(m->vbo.size()+instances.size(), sizeof(glm::mat4)/sizeof(GLfloat), GL_FLOAT, GL_FALSE, 0, (void*)0);
  glVertexAttribDivisor(m->vbo.size()+instances.size(), 1);
  instances.push_back(instance);*/
}


void Instance::updateBuffer(std::vector<glm::mat4>& buf, int index){
  glBindVertexArray(m->vao);
  glBindBuffer(GL_ARRAY_BUFFER, instances[index]);  //Bind Instance Buffer and Data
  if(buf.size() != SIZE)  glBufferData(GL_ARRAY_BUFFER, buf.size()*sizeof(glm::mat4), &buf[0], GL_STATIC_DRAW);
  else                    glBufferSubData(GL_ARRAY_BUFFER, 0, SIZE*sizeof(glm::mat4), &buf[0]);
  SIZE = buf.size();
}

void Instance::render(GLenum mode){
  glBindVertexArray(m->vao);
  Model *model = dynamic_cast<Model *>(m);
  if(model && model->indexed)
  {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model->ibo);
    glDrawElementsInstanced(mode, model->SIZE, GL_UNSIGNED_INT, 0, SIZE);
  }
  else
    glDrawArraysInstanced(mode, 0, m->SIZE, SIZE); //Instanced render
}
