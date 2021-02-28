#pragma once
#include "field_2d.h"
#include "terrain.h"
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