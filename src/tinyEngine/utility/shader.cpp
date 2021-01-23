#include "shader.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <iosfwd>
#include <string>
#include <glm/glm.hpp>
#include <sstream>
#include "../utility.h"
void Shader::setup(slist _s){
  boost::filesystem::path data_dir("");
  std::vector<std::string> s = _s;
  for (std::string &path : s)
  {
    path = base_shader_path + path;
  }
  if(s.size() == 2){
    vertexShader   = addProgram((data_dir/s[0]).string(), GL_VERTEX_SHADER);
    fragmentShader = addProgram((data_dir/s[1]).string(), GL_FRAGMENT_SHADER);
  }
  else if(s.size() == 3){
    vertexShader   = addProgram((data_dir/s[0]).string(), GL_VERTEX_SHADER);
    geometryShader = addProgram((data_dir/s[1]).string(), GL_GEOMETRY_SHADER);
    fragmentShader = addProgram((data_dir/s[2]).string(), GL_FRAGMENT_SHADER);
  }
  else std::cout<<"Number of shaders not recognized."<<std::endl;
}

int Shader::addProgram(std::string fileName, GLenum shaderType){
  char* src; int32_t size;
  std::string result = readGLSLFile(fileName, size);
  src = const_cast<char*>(result.c_str());

  int shaderID = glCreateShader(shaderType);
  glShaderSource(shaderID, 1, &src, &size);
  compile(shaderID);

  return shaderID;
}

void Shader::compile(GLuint shader){
  glCompileShader(shader);
  int success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if(success) glAttachShader(program, shader);
  else        error(shader, true);
}

void Shader::link(){
  glLinkProgram(program);
  int success;
  glGetProgramiv(program, GL_LINK_STATUS, &success);
  if(!success) error(program, false);
}

void Shader::error(GLuint s, bool t){
  int m;
  if(t) glGetShaderiv(s, GL_INFO_LOG_LENGTH, &m);
  else glGetProgramiv(s, GL_INFO_LOG_LENGTH, &m);
  char* l = new char[m];
  if(t) glGetShaderInfoLog(s, m, &m, l);
  else glGetProgramInfoLog(s, m, &m, l);
  std::cout<<"Linker Error: "<<l<<std::endl;
  delete[] l;
}

void Shader::use(){
  boundtextures = 0;
  glUseProgram(program);
}

std::string Shader::readGLSLFile(std::string file, int32_t &size){
  std::ifstream t;
  std::string fileContent;

  t.open(file);     //Read GLSL File Contents
  if(t.is_open()){
    std::stringstream buffer;
    buffer << t.rdbuf();
    fileContent = buffer.str();
  }
  else std::cout<<"File "<<file<<" opening failed"<<std::endl;
  t.close();

  size = fileContent.length();  //Set the Size
  return fileContent;
}

/* Shader Storage Buffer Objects */

void Shader::addBuffer(std::string name){
  unsigned int b;
  glGenBuffers(1, &b);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, b);
  glShaderStorageBlockBinding(program, glGetProgramResourceIndex(program, GL_SHADER_STORAGE_BLOCK, name.c_str()), ssbo.size());
  glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ssbo.size(), b);
  sbpi[name] = ssbo.size();
  ssbo[name] = b;
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}
void Shader::texture(std::string name, const Texture& t)
{
  glActiveTexture(GL_TEXTURE0 + boundtextures);
  glBindTexture(t.type, t.texture);
  uniform(name, boundtextures++);
}

