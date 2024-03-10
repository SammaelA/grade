#pragma once
#include "common_utils/field_2d.h"
#include "graphics_utils/terrain.h"
#include "save_utils/serialization.h"

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
    float2 pos;
    int roots_count;
};
class PlanarShadowsMap : public Field_2d
{
public:
    PlanarShadowsMap(float3 pos, float2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void set_occluder(float3 position, float base_val, float r, float pow);
    void clear();
    void add(PlanarShadowsMap &src);
    void set(PlanarShadowsMap &src);
    void add_body(Body *b, float opacity = 1e9, bool solid = true);
};
class GroveMask : public Field_2d
{
public:
    friend class boost::serialization::access;

    GroveMask(float3 pos, float2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    GroveMask() : Field_2d() {};
    void set_round(float r);
    void set_round_min(float r, float val);
    void set_square(float x, float z);
private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<Field_2d>(*this);
    }
};
class HabitabilityMap : public Field_2d
{
public:
    HabitabilityMap(float3 pos, float2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void create(Heightmap &h, GroveMask &mask);
};
class DensityMap : public Field_2d
{
public:
    DensityMap(float3 pos, float2 size, float cell_size) : Field_2d(pos, size, cell_size) {};
    void clear();
    void create(HabitabilityMap &hm, PlanarShadowsMap &psm);
    void choose_places_for_seeds(int count, std::vector<Seed> &seeds);
private:
    double calc_sum();
};