#pragma once
#include "common_utils/parameter.h"
#include "core/tree.h"
#include "graphics_utils/terrain.h"
#include "generation/generation_settings.h"
#include <atomic>

struct LightVoxelsCube;
class AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) = 0;

    virtual bool iterate(LightVoxelsCube &voxels) { return false;};//return true if everything is finished
    virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) {};
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) {};
    virtual bool iteration_method_implemented() {return false;}
    virtual bool use_voxels_for_generation() {return false;}
    virtual void set_seed(int seed) {};
    static void set_joints_limit(int lim) {joints_limit = lim;}
protected:
    static int joints_limit;
    static std::atomic<int> branch_next_id, tree_next_id; 
};