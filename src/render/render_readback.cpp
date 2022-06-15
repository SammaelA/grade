#include "render_readback.h"
#include "cities_generator/global.h"
#include "graphics_utils/texture_manager.h"

RenderReadback::RenderReadback():
screen_to_world_pos({"screen_to_world_pos.comp"},{})
{
    glGenBuffers(1, &results_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 13, results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec4), NULL, GL_STREAM_READ);
}
RenderReadback::~RenderReadback()
{
    glDeleteBuffers(1, &results_buf);
}
glm::vec4 RenderReadback::get_world_pos(glm::vec2 screen_pos, GLuint world_pos_tex, GLuint color_tex)
{
    screen_to_world_pos.use();
    screen_to_world_pos.texture("world_pos_tex", world_pos_tex);
    screen_to_world_pos.texture("color_tex", color_tex);
    screen_to_world_pos.uniform("screen_pos", screen_pos);

    glDispatchCompute(1, 1, 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_WRITE);
    glm::vec4 *vec_ptr = (glm::vec4*)(ptr);
    glm::vec4 res = vec_ptr[0];
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    return res;
}