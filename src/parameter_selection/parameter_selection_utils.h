#pragma once
#include "core/scene.h"
#include "graphics_utils/volumetric_occlusion.h"

namespace ps_utils
{
    void gen_tree(LightVoxelsCube &voxels, TreeTypeData *type, Tree *tree);
    void gen_tree_task(int start_n, int stop_n, LightVoxelsCube *vox, std::vector<TreeTypeData> *types, Tree *trees);
};