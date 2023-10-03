#include "resources.h"
#include "common_utils/utility.h"

int buffers_cnt = 0;
int framebuffers_cnt = 0;
int vertex_arrays_cnt = 0;
int queries_cnt = 0;

unsigned int create_buffer()
{
  buffers_cnt++;
  GLuint id = 0;
  glGenBuffers(1, &id);
  return id;
}
void delete_buffer(unsigned int id)
{
  if (id != ~0u && id != 0u)
  {
    buffers_cnt--;
    glDeleteBuffers(1, &id);
  }
}

unsigned int create_framebuffer()
{
  framebuffers_cnt++;
  GLuint id = 0;
  glGenFramebuffers(1, &id);
  return id;
}
void delete_framebuffer(unsigned int id)
{
  if (id != ~0u && id != 0u)
  {
    framebuffers_cnt--;
    glDeleteFramebuffers(1, &id);
  }
}

unsigned int create_vertex_array()
{
  vertex_arrays_cnt++;
  GLuint id = 0;
  glGenVertexArrays(1, &id);
  return id;
}
void delete_vertex_array(unsigned int id)
{
  if (id != ~0u && id != 0u)
  {
    vertex_arrays_cnt--;
    glDeleteVertexArrays(1, &id);
  }
}

unsigned int create_query()
{
  queries_cnt++;
  GLuint id = 0;
  glGenQueries(1, &id);
  return id;
}
void delete_query(unsigned int id)
{
  if (id != ~0u && id != 0u)
  {
    queries_cnt--;
    glDeleteQueries(1, &id);
  }
}

void print_FB_status(GLuint status)
{
  if (status == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
    debugl(9, "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT\n");
  else if (status == GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS)
    debugl(9, "GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS\n");
  else if (status == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT)
    debugl(9, "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT\n");
  else if (status == GL_FRAMEBUFFER_UNSUPPORTED)
    debugl(9, "GL_FRAMEBUFFER_UNSUPPORTED\n");
  else
    debugl(9, "GL_FRAMEBUFFER_INCOMPLETE %#010x\n", status);
}