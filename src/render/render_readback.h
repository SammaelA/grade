#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/shader.h"
#include "render_readback_data.h"

class RenderReadback
{
public:
    RenderReadback();
    ~RenderReadback();
    float4 get_world_pos(float2 screen_pos, Texture &world_pos_tex, Texture &color_tex);//w component is < 0 if there is nothing at screen_pos
private:
    GLuint results_buf;
    Shader screen_to_world_pos;
};