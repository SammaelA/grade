#pragma once
#include <GL/glew.h>
#include <string>
#include <initializer_list>
#include <unordered_map>
#include <vector>
#include <iostream>
#include "glm/glm.hpp"
#include "texture.h"
class Shader{
using slist = std::initializer_list<std::string>;
public:
  const std::string base_shader_path = "/home/sammael/study/bit_bucket/grade/shaders/";
  template<typename... Args>
  Shader(slist shaders, slist in){
    program = glCreateProgram();        //Generate Shader
    setup(shaders);                     //Add Individual Shaders
    for(auto &n : in)                   //Add all Attributes of Shader
      glBindAttribLocation(program, &n - in.begin(), n.c_str());
    link();                             //Link the shader program!
  }

  Shader(slist shaders, slist in){
    program = glCreateProgram();        //Generate Shader
    setup(shaders);                     //Add Individual Shaders
    for(auto &n : in)                   //Add all Attributes of Shader
      glBindAttribLocation(program, &n - in.begin(), n.c_str());
    link();                             //Link the shader program!
  }

  Shader(slist shaders, slist in, slist buf):Shader(shaders, in){
    for(auto&b : buf) addBuffer(b);     //Add Storage Buffers to Shader
  }

  ~Shader(){
    glDeleteProgram(program);
    glDeleteShader(fragmentShader);
    glDeleteShader(geometryShader);
    glDeleteShader(vertexShader);
    glDeleteShader(computeShader);
  }

  GLuint program;   //Shader Program ID
  GLuint vertexShader, geometryShader, fragmentShader, computeShader;
  int boundtextures;

  std::unordered_map<std::string, GLuint> ssbo; //SSBO Storage
  std::unordered_map<std::string, GLuint> sbpi; //Shader Binding Point Index

  void setup(slist shaders);
  int  addProgram(std::string fileName, GLenum shaderType);  //General Shader Addition
  std::string readGLSLFile(std::string fileName, int32_t &size); //Read File
  void compile(GLuint shader);  //Compile and Add File
  void link();                  //Link the entire program
  void error(GLuint s, bool t); //Get Compile/Link Error
  void use();                   //Use the program
  void addBuffer(std::string name);                          //General Buffer Addition

  void texture(std::string name, const Texture& t);
  void texture(std::string name, GLuint t);

/* Uniform Setters */

  void  uniform(std::string name, const bool u){
  glUniform1i(glGetUniformLocation(program, name.c_str()), u); }

  void  uniform(std::string name, const int u){
  glUniform1i(glGetUniformLocation(program, name.c_str()), u); }
  
  void  uniform(std::string name, const uint u){
  glUniform1ui(glGetUniformLocation(program, name.c_str()), u); }
  
  void  uniform(std::string name, const float u){
  glUniform1f(glGetUniformLocation(program, name.c_str()), u); }

  void  uniform(std::string name, const double u){ //GLSL Intrinsically Single Precision
  glUniform1f(glGetUniformLocation(program, name.c_str()), (float)u); }

  void  uniform(std::string name, const glm::vec2 u){
  glUniform2fv(glGetUniformLocation(program, name.c_str()), 1, &u[0]); }

  void  uniform(std::string name, const glm::vec3 u){
  glUniform3fv(glGetUniformLocation(program, name.c_str()), 1, &u[0]); }

  void  uniform(std::string name, const float (&u)[3]){
  glUniform3fv(glGetUniformLocation(program, name.c_str()), 1, &u[0]); }

  void  uniform(std::string name, const float (&u)[4]){
  glUniform4fv(glGetUniformLocation(program, name.c_str()), 1, &u[0]); }

  void  uniform(std::string name, const glm::vec4 u){
  glUniform4fv(glGetUniformLocation(program, name.c_str()), 1, &u[0]); }

  void  uniform(std::string name, const glm::mat3 u){
  glUniformMatrix3fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, &u[0][0]); }

  void  uniform(std::string name, const glm::mat4 u){
  glUniformMatrix4fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, &u[0][0]); }

  void  uniform(std::string name, const std::vector<glm::mat4> u){
  glUniformMatrix4fv(glGetUniformLocation(program, name.c_str()), u.size(), GL_FALSE, &u[0][0][0]); }
};