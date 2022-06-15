#pragma once
#include "tinyEngine/texture.h"
#include "tinyEngine/shader.h"
#include "render_readback_data.h"

class RenderReadback
{
public:
    RenderReadback();
    ~RenderReadback();
    glm::vec4 get_world_pos(glm::vec2 screen_pos, GLuint world_pos_tex, GLuint color_tex);//w component is < 0 if there is nothing at screen_pos
private:
    GLuint results_buf;
    Shader screen_to_world_pos;
};