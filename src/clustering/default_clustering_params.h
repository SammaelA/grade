#pragma once
#include "save_utils/blk.h"

enum ClusteringStep
{
    TRUNKS,
    BRANCHES,
    TREES,
    CLUSTERING_STEPS_COUNT
};

extern Block trunks_default_params;
extern Block branches_default_params;
extern Block trees_default_params;
extern ClusteringStep current_clustering_step;
Block &get_default_block();
void load_default_blocks();