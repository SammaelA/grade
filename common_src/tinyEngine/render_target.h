#pragma once
#include "texture.h"

class RenderTarget
{
public:
    RenderTarget(){}
    ~RenderTarget();
    bool create(int w, int h);
    void target();

    void set_clear_color(float4 cl_color) {clear_color = cl_color;}
    Texture get_tex(){ return tex;}
private:
    float4 clear_color = float4(0,0,0,1);
    int width, height;
    GLuint frBuffer;
    Texture tex;
    GLuint texFmt = GL_RGBA8;
};