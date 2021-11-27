#pragma once
#include "volumetric_occlusion.h"

struct DebugTransferData
{
    std::vector<LightVoxelsCube *> debug_voxels;
};
struct DebugTransferSettings
{
    int save_detailed_voxels_count = 0;
    int save_small_voxels_count  = 0;
};

extern DebugTransferData debugTransferData;
extern DebugTransferSettings debugTransferSettings;