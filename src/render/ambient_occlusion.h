#pragma once
#include "tinyEngine/shader.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/engine.h"

class HBAORenderer
{
public:
    HBAORenderer();
    ~HBAORenderer();
    void render(float fovRad, float width, float height, GLuint viewPosTex);
    void create(int w, int h);
    Texture &get_tex() {return aoTex;}
private:
    Texture noise;
    PostFx shader;
    float4 clear_color = float4(0,0,0,0);
    int width, height;
    GLuint frBuffer;
    Texture aoTex;
    GLuint aoTexFmt = GL_R16F; 
};
