#include "../camera.h"
#include "utility/texture.h"
#include "utility/shader.h"
#include "utility/postfx.h"
#pragma once
class ShadowMap
{
    public:
    void use(DirectedLight &light);
    void blur();
    glm::mat4 &get_transform();
    glm::mat4 &get_view(){return view;}
    glm::mat4 &get_projection(){return projection;}
    GLuint getTex() {return VSMdepthTex;}
    void create(int w, int h);
    ~ShadowMap();

    GLuint depthMapFBO, depthMap;
    GLuint VSMdepthTex, VSMdepthTexTemp;
    int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;
    glm::mat4 view;
    glm::mat4 projection;
    glm::mat4 viewproj;
    Camera shadow_camera;
    PostFx *postFx = nullptr;
};