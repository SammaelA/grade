#pragma once
#include "tinyEngine/camera.h"
#include "tinyEngine/texture.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/postfx.h"

class ShadowMap
{
    public:
    void use(DirectedLight &light);
    void start_trans_pass();
    void finish_trans_pass();
    void blur();
    float4x4 &get_transform();
    float4x4 &get_view(){return view;}
    float4x4 &get_projection(){return projection;}
    Texture &getTex() {return VSMdepthTex;}
    void create(int w, int h);
    ~ShadowMap();

    GLuint depthMapFBO;
    Texture depthMap, VSMdepthTex, VSMdepthTexTemp, srcDepthTex;
    int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;
    float4x4 view;
    float4x4 projection;
    float4x4 viewproj;
    Camera shadow_camera;
    PostFx *postFx = nullptr;
};