#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include "tinyEngine/utility/texture.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/save_utils/saver.h"
class BillboardCloudRenderer;
class TextureAtlas
{
public:
    friend class BillboardCloudRaw;
    friend bool saver::save(FILE *f, TextureAtlas &t); 
    friend bool saver::load(FILE *f, TextureAtlas &t); 
    TextureAtlas();
    TextureAtlas(const TextureAtlas &a);
    TextureAtlas(int w, int h, int layers = 2);
    ~TextureAtlas();
    void set_clear_color(glm::vec4 color)
    {
        clearColor = color;
    }
    TextureAtlas &operator=(TextureAtlas &atlas);
    void set_grid(int w, int h);
    int add_tex();
    void process_tc(int num, glm::vec3 &tc);
    bool target(int num);
    glm::mat4 tex_transform(int num);
    glm::ivec4 get_sizes() { return glm::ivec4(width, height, gridWN, gridHN); }
    bool clear();
    Texture &tex() { return colorTex; }
    void gen_mipmaps();
    int layers_count() {return layers;}
private:
    bool bind(int layer);
    int curNum = 0;
    int width = 0;
    int height = 0;
    int layers = 0;
    int gridWN = 0;
    int gridHN = 0;
    bool isGrid = false;
    glm::vec4 clearColor;
    GLuint fbo;
    Texture colorTex;
    Shader mipMapRenderer;
    Shader copy;
};