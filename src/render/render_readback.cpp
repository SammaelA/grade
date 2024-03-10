#include "render_readback.h"
#include "tinyEngine/engine.h"

RenderReadback::RenderReadback():
screen_to_world_pos({"screen_to_world_pos.comp"},{})
{
    results_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 13, results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float4), NULL, GL_STREAM_READ);
}
RenderReadback::~RenderReadback()
{
    delete_buffer(results_buf);
}
float4 RenderReadback::get_world_pos(float2 screen_pos, Texture &world_pos_tex, Texture &color_tex)
{
    screen_to_world_pos.use();
    screen_to_world_pos.texture("world_pos_tex", world_pos_tex);
    screen_to_world_pos.texture("color_tex", color_tex);
    screen_to_world_pos.uniform("screen_pos", screen_pos);

    glDispatchCompute(1, 1, 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_WRITE);
    float4 *vec_ptr = (float4*)(ptr);
    float4 res = vec_ptr[0];
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    return res;
}