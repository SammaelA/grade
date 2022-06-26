#pragma once
#include "shader.h"
#include "texture.h"
class DefferedTarget
{
public:
    DefferedTarget(){}
    ~DefferedTarget();
    bool create(int w, int h);
    void target();

    void set_clear_color(glm::vec4 cl_color) {clear_color = cl_color;}
    Texture &get_color(){ return colorTex;}
    Texture &get_normals(){ return normalsTex;}
    Texture &get_view_pos(){ return viewPosTex;}
    Texture &get_world_pos(){ return worldPosTex;}
    glm::vec2 size(){ return glm::vec2(width,height);}
private:
    glm::vec4 clear_color;
    int width, height;
    GLuint frBuffer;
    Texture colorTex, normalsTex, viewPosTex, depthTex, worldPosTex;
    GLuint colorTexFmt = GL_RGBA8, normalsTexFmt = GL_RGBA16F, viewPosTexFmt = GL_RGB16F, worldPosTexFmt = GL_RGB16F;
};