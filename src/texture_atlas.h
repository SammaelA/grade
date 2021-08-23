#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include "tinyEngine/utility/texture.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/bit_vector.h"
#include "tinyEngine/save_utils/saver.h"
class BillboardCloudRenderer;
class TextureAtlas;
enum Channel
{
    R,
    G,
    B,
    A
};
class TextureAtlasRawData
{
public:
    unsigned char get_pixel_uc(int w, int h, Channel chan, int tex_id);
    unsigned char get_pixel_uc_safe(int w, int h, Channel chan, int tex_id);
    float get_pixel(int w, int h, Channel chan, int tex_id);
    TextureAtlasRawData();
    TextureAtlasRawData(TextureAtlas &atlas);
    ~TextureAtlasRawData();
    void clear();
    
    int get_w() { return w; }
    int get_h() { return h; }
    int get_slices() { return slices; }
    int get_layers() { return layers; }
    
    bool is_valid() { return valid; }
private:
    bool valid = false;
    int w = 0,h = 0, layers = 0, slices = 0;
    int gridWN = 0;
    int gridHN = 0;
    int slice_h;
    int slice_w;
    unsigned char *raw_data = nullptr;
};
class TextureAtlas : public Countable
{
public:
    friend class BillboardCloudRaw;
    friend class TextureAtlasRawData;
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
    void set_grid(int w, int h, bool resizable = true);
    int add_tex();
    void remove_tex(int pos);
    void process_tc(int num, glm::vec3 &tc);
    glm::vec4 tc_transform(int num);
    bool target(int num, int type);
    glm::mat4 tex_transform(int num);
    glm::ivec4 get_sizes() { return glm::ivec4(width, height, gridWN, gridHN); }
    bool clear();
    Texture &tex(int type);
    void gen_mipmaps();
    int layers_count() {return layers;}
    int tex_count() {return 2;}
    int capacity() { return gridWN*gridHN*layers;}
    bool is_valid() {return valid && (capacity() > 0);}
    void destroy();
private:
    bool bind(int layer, int type);
    int new_layers_count();
    void increase_capacity_tex(Texture &t, int layers);
    void increase_capacity();
    int curNum = 0;
    int width = 0;
    int height = 0;
    int layers = 0;
    int gridWN = 0;
    int gridHN = 0;
    bool isGrid = false;
    bool resizable = false;
    bool valid = false;
    glm::vec4 clearColor;
    GLuint fbo;
    Texture colorTex;
    Texture normalTex;
    Shader mipMapRenderer;
    Shader copy;
    BitVector occupied;
};