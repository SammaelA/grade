#ifndef SHADER_H
#define SHADER_H

#include <string>

#include <GL/glew.h> 

class Shader
{
public:
    
    GLuint Program;
    const std::string name;
    
    Shader(const GLchar* vertexPath, const GLchar* fragmentPath, std::string shaderName);
    
    void Use();
};

#endif