#pragma once
#include "field_2d.h"
#include "terrain.h"
struct GroveGenerationData;
struct Tree;
struct Branch;
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
};
class GroveMask : public Field_2d
{
public:
    GroveMask(glm::vec3 pos, glm::vec2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void set_round(float r);
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
class Seeder
{
public:
    Seeder(GroveGenerationData &ggd, float cell_size, Heightmap *h);
    void recalcuate_shadows(Tree *trees, int count);
    void add_tree_shadow(Tree &t);
    void choose_places_for_seeds(int count, std::vector<Seed> &seeds);
    Heightmap *heightmap;
private:
    void recalculate_planar_shadows(Branch *b, PlanarShadowsMap &psm, int level);
    int joints_count(Branch *b);
    GroveMask mask;
    HabitabilityMap hm;
    PlanarShadowsMap psm;
    DensityMap dsm;
};