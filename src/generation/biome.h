#pragma once
#include <vector>
#include "save_utils/blk.h"
#include "common_utils/utility.h"
#include "common_utils/bbox.h"

struct Biome
{
    struct PatchDesc
    {
        float size = 25;
        float size_std_dev = 10;
        float density = 1;
        float density_std_dev = 0.5;
        float coverage_part = 0.0;
        float push = 1;

        std::vector<std::pair<int, float>> types;
    };
    struct VegClass
    {
        float main_density = -1;
        std::vector<std::pair<int, float>> main_types;
        std::vector<PatchDesc> patch_descs;
    };
    
    int id = -1;
    VegClass trees;
    VegClass plants;
    VegClass grass;

    void load_from_blk(Block &b);
};

class BiomeMap
{
public:
    typedef unsigned char biome_type_t; 

    BiomeMap();
    ~BiomeMap();
    void create(AABB2D bbox, float pixel_size);
    int pixels_w() {return w;}
    int pixels_h() {return h;}
    AABB2D borders() {return bbox;}
    int get(int w, int h);
    int get(glm::vec2 pos);
    int get(glm::vec3 pos);
    
    void set_rect(AABB2D box, int id);
    void set_round(glm::vec2 pos, float inner_r, float outer_r, int id);
    void save_as_image(std::string name = "biome_map_debug");
private:

   biome_type_t *data = nullptr;
   int w,h;
   AABB2D bbox;
   float pixel_size;
};