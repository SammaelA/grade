#pragma once
#include "tinyEngine/resources.h"
#include <glm/glm.hpp>
#include <vector>
#include "tinyEngine/texture.h"
#include "tinyEngine/shader.h"
#include "common_utils/bit_vector.h"
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
    unsigned char *get_raw_data() { return raw_data;}
    int get_slice_size(int tex_id);
    void get_slice(int tex_id, unsigned char *data, int *slice_w = nullptr, int *slice_h = nullptr);
    float get_pixel(int w, int h, Channel chan, int tex_id);
    TextureAtlasRawData();
    TextureAtlasRawData(const TextureAtlas &atlas);
    ~TextureAtlasRawData();
    void clear();
    
    int get_w() { return w; }
    int get_h() { return h; }
    int get_slices() { return slices; }
    int get_layers() { return layers; }
    
    bool is_valid() { return valid; }
private:
    int get_pixel_pos(int w, int h, int tex_id);
    bool valid = false;
    int w = 0,h = 0, layers = 0, slices = 0;
    int gridWN = 0;
    int gridHN = 0;
    int slice_h = 0;
    int slice_w = 0;
    unsigned char *raw_data = nullptr;
};
class TextureAtlas : public Countable
{
public:
    friend class BillboardCloudRaw;
    friend class TextureAtlasRawData;
    TextureAtlas();
    TextureAtlas(const TextureAtlas &a);
    TextureAtlas(int w, int h, int layers, int mip_levels = 5, bool need_normals = false);
    ~TextureAtlas();
    void set_clear_color(glm::vec4 color)
    {
        clearColor = color;
    }
    TextureAtlas &operator=(const TextureAtlas &atlas);
    void set_grid(int w, int h, bool resizable = true);
    int add_tex();
    void remove_tex(int pos);
    void process_tc(int num, glm::vec3 &tc) const;
    glm::vec4 tc_transform(int num) const;
    void pixel_offsets(int num, glm::ivec3 &tc) const;
    void pixel_offsets(int num, glm::uvec4 &tc) const;
    bool target(int num, int type);
    bool target_slice(int num, int type);
    glm::mat4 tex_transform(int num) const;
    glm::ivec4 get_sizes() const { return glm::ivec4(width, height, gridWN, gridHN); }
    glm::ivec2 get_slice_size() const {return glm::ivec2(width/gridWN, height/gridHN); }
    bool clear();
    Texture &tex(int type);
    const Texture &color() const {return colorTex;}
    const Texture &normal() const {return normalTex;}
    void gen_mipmaps(std::string mipmap_shader_name = "mipmap_render.fs");
    int layers_count() const {return layers;}
    int tex_count() const {return 2;}
    int capacity() const { return gridWN*gridHN*layers;}
    bool is_valid() const {return valid && (capacity() > 0);}
    std::vector<int> get_all_valid_slices_ids();
    void destroy();

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
    GLuint fbo = 0;
    Texture colorTex;
    Texture normalTex;
    Shader mipMapRenderer;
    Shader copy;
    BitVector occupied;
};