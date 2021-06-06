#pragma once
#include "shader.h"
#include "target.h"

class DefferedTarget
{
public:
    DefferedTarget(){}
    ~DefferedTarget();
    bool create(int w, int h);
    void target();

    void set_clear_color(glm::vec4 cl_color) {clear_color = cl_color;}
    GLuint get_color(){ return colorTex;}
    GLuint get_normals(){ return normalsTex;}
    GLuint get_view_pos(){ return viewPosTex;}
    GLuint get_world_pos(){ return worldPosTex;}
private:
    void print_FB_status(GLuint status);
    glm::vec4 clear_color;
    int width, height;
    GLuint frBuffer;
    GLuint colorTex, normalsTex, viewPosTex, depthTex, worldPosTex;
    GLuint colorTexFmt = GL_RGBA8, normalsTexFmt = GL_RGBA16F, viewPosTexFmt = GL_RGB16F, worldPosTexFmt = GL_RGB16F;
};