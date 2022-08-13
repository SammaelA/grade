#pragma once
#include "common_utils/parameter.h"
#include "core/tree.h"
#include "graphics_utils/terrain.h"
#include "generation/generation_task.h"
#include <atomic>

struct LightVoxelsCube;
class AbstractTreeGenerator
{
public:
    virtual bool iterate(LightVoxelsCube &voxels) { return false;};//return true if everything is finished
    virtual void plant_tree(glm::vec3 pos, const TreeTypeData *type) = 0;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) = 0;
    virtual bool use_voxels_for_generation() {return false;}
    virtual void set_seed(int seed) {};
    static void set_joints_limit(int lim) {joints_limit = lim;}
protected:
    static int joints_limit;
    static std::atomic<int> branch_next_id, tree_next_id; 
};