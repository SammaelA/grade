#pragma once
#include "tinyEngine/shader.h"
#include "tinyEngine/postfx.h"
#include "tinyEngine/camera.h"
#include "app.h"
#include "graphics_utils/texture_manager.h"

class HBAORenderer
{
public:
    HBAORenderer();
    ~HBAORenderer();
    void render(AppContext &ctx, GLuint viewPosTex);
    void create(int w, int h);
    GLuint get_tex() {return aoTex;}
private:
    Texture noise;
    PostFx shader;
    glm::vec4 clear_color = glm::vec4(0,0,0,0);
    int width, height;
    GLuint frBuffer;
    GLuint aoTex;
    GLuint aoTexFmt = GL_R16F; 
};
