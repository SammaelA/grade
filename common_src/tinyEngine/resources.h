#pragma once
#include <GL/glew.h> 
#include <string>

unsigned int create_buffer();
void delete_buffer(unsigned int id);

unsigned int create_framebuffer();
void delete_framebuffer(unsigned int id);

unsigned int create_vertex_array();
void delete_vertex_array(unsigned int id);

unsigned int create_query();
void delete_query(unsigned int id);

void print_FB_status(GLuint status);

void checkGLErrors(const std::string &text, bool clean = false);