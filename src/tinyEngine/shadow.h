#include "../camera.h"
#include "utility/texture.h"
#include "utility/shader.h"
#pragma once
class ShadowMap
{
    public:
    void use(DirectedLight &light);
    glm::mat4 get_transform();
    GLuint getTex() {return depthMap;}
    void create(int w, int h);
    ~ShadowMap() {glDeleteTextures(1,&depthMap);}

    GLuint depthMapFBO, depthMap;
    Texture *tm;
    int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;
    glm::mat4 viewproj;
    Camera shadow_camera;
};