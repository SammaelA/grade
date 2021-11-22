#pragma once
#include "core/grass.h"
#include "common_utils/bbox.h"
#include "graphics_utils/terrain.h"
#include <map>
#include <vector>
#include <list>

struct Cell;
struct Block;
struct LightVoxelsCube;
class GrassGenerator
{
public:
    void prepare_grass_patches(std::vector<Cell> &cells, int cells_x, int cells_y);
    void set_grass_types(const std::map<std::string, GrassType> &grass_types, Block &grass_settings);
    void generate_grass_in_cell(Cell &cell, LightVoxelsCube *occlusion);
    void pack_all_grass(GrassPacked &grass_packed, Heightmap &h);
private:
    struct GrassPatch
    {
        int grass_type;//position in used_grass_types vector
        Sphere2D sphere;
        int instances_max_cnt;
        std::list<Sphere2D> grass_instances;
    };
    std::vector<GrassPatch> grass_patches;
    std::vector<GrassType> used_grass_types;
    std::vector<float> grass_quantity;//for every type above

    void add_patch(std::vector<Cell> &cells, int cells_x, int cells_y, int cell_id, int type_id);
};