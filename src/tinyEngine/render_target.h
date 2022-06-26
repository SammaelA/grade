#pragma once
#include "texture.h"

class RenderTarget
{
public:
    RenderTarget(){}
    ~RenderTarget();
    bool create(int w, int h);
    void target();

    void set_clear_color(glm::vec4 cl_color) {clear_color = cl_color;}
    Texture get_tex(){ return tex;}
private:
    glm::vec4 clear_color;
    int width, height;
    GLuint frBuffer;
    Texture tex;
    GLuint texFmt = GL_RGBA8;
};