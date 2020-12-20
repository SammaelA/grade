#pragma once
#include <GL/glew.h>  
#include <glm/glm.hpp>
#include <vector>
#include "tinyEngine/utility/texture.h"
#include "tinyEngine/utility/shader.h"

class TextureAtlas
{
    public:
    TextureAtlas(int w, int h);
    ~TextureAtlas()
    {
        glDeleteFramebuffers(1, &fbo);
    }
    void set_clear_color(glm::vec4 color)
    {
        clearColor = color;
    }
    void set_grid(int w, int h);
    int add_tex();
    void process_tc(int num, glm::vec2 &tc);
    void process_tc_array(int num, std::vector<float> &tc);
    bool target(int num);
    glm::mat4 tex_transform(int num);
    glm::ivec4 get_sizes() {return glm::ivec4(width,height,gridWN, gridHN);}
    bool clear();
    Texture &tex() {return colorTex;}
    void gen_mipmaps();
    private:
    bool bind();
    int curNum = 0;
    int width = 0;
    int height = 0;
    int gridWN = 0;
    int gridHN = 0;
    bool isGrid = false;
    glm::vec4 clearColor;
    GLuint fbo;  
    Texture colorTex;
    Shader mipMapRenderer;
    Shader copy;
};