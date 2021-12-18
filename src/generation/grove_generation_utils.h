#pragma once
#include "common_utils/field_2d.h"
#include "graphics_utils/terrain.h"

struct GroveGenerationData;
struct Tree;
struct Branch;
class Body;
namespace mygen
{
class Tree;
class Branch;
}
struct Seed
{
    glm::vec2 pos;
    int roots_count;
};
class PlanarShadowsMap : public Field_2d
{
public:
    PlanarShadowsMap(glm::vec3 pos, glm::vec2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void set_occluder(glm::vec3 position, float base_val, float r, float pow);
    void clear();
    void add(PlanarShadowsMap &src);
    void set(PlanarShadowsMap &src);
    void add_body(Body *b, float opacity = 1e9, bool solid = true);
};
class GroveMask : public Field_2d
{
public:
    GroveMask(glm::vec3 pos, glm::vec2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void set_round(float r);
    void set_round_min(float r, float val);
    void set_square(float x, float z);
};
class HabitabilityMap : public Field_2d
{
public:
    HabitabilityMap(glm::vec3 pos, glm::vec2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void create(Heightmap &h, GroveMask &mask);
};
class DensityMap : public Field_2d
{
public:
    DensityMap(glm::vec3 pos, glm::vec2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void clear();
    void create(HabitabilityMap &hm, PlanarShadowsMap &psm);
    void choose_places_for_seeds(int count, std::vector<Seed> &seeds);
private:
    double calc_sum();
};
class Seeder : Countable
{
public:
    Seeder(GroveGenerationData &ggd, float cell_size, const Heightmap *h);
    Seeder(glm::vec3 pos, glm::vec3 size, float cell_size, const Heightmap *h);
    void recalcuate_shadows(Tree *trees, int count);
    void add_tree_shadow(Tree &t);

    void recalcuate_shadows(mygen::Tree *trees, int count);
    void add_tree_shadow(mygen::Tree &t);

    void choose_places_for_seeds(int count, std::vector<Seed> &seeds);
    void add_body(Body *b, float opacity = 1e9, bool solid = true);
    Heightmap *heightmap;
private:
    void recalculate_planar_shadows(Branch *b, PlanarShadowsMap &psm, int level);
    int joints_count(Branch *b);

    void recalculate_planar_shadows(mygen::Branch *b, PlanarShadowsMap &psm, int level);
    int joints_count(mygen::Branch *b);

    GroveMask mask;
    HabitabilityMap hm;
    PlanarShadowsMap psm;
    PlanarShadowsMap const_psm;
    DensityMap dsm;
};