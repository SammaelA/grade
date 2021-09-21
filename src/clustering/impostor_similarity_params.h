#pragma once
#include "../tinyEngine/save_utils/blk.h"
#include "default_clustering_params.h"

struct ImpostorSimilarityParams
{
    int impostor_similarity_slices = 4;
    int impostor_metric_level_from = 0;
    int impostor_metric_level_to = 1000;
    float leaf_size_mult = 1;
    float wood_size_mult = 1;
    void load_from_block(Block *b)
    {
        impostor_similarity_slices = b->get_int("impostor_similarity_slices",impostor_similarity_slices);
        impostor_metric_level_from = b->get_int("impostor_metric_level_from",impostor_metric_level_from);
        impostor_metric_level_to = b->get_int("impostor_metric_level_to",impostor_metric_level_to);
        leaf_size_mult = b->get_double("leaf_size_mult",leaf_size_mult);
        wood_size_mult = b->get_double("wood_size_mult",wood_size_mult);
    }
    void load(Block *b)
    {
        if (!b)
            return;
        
        Block &def = get_default_block();
        load_from_block(&def);
        load_from_block(b);
    }
};