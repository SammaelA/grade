#pragma once
#include "common_utils/blk.h"
#include "default_clustering_params.h"

struct ImpostorSimilarityParams
{
    int impostor_similarity_slices = 4;
    int impostor_metric_level_from = 0;
    int impostor_metric_level_to = 1000;
    int impostor_texture_size = Quality::LOW_AS_F;
    float leaf_size_mult = 1;
    float wood_size_mult = 1;
    float size_diff_factor = 0.75;
    float size_diff_tolerance = 0.1;
    float leaves_opacity = 0.4;
    void load_from_block(Block *b)
    {
        impostor_similarity_slices = b->get_int("impostor_similarity_slices",impostor_similarity_slices);
        impostor_metric_level_from = b->get_int("impostor_metric_level_from",impostor_metric_level_from);
        impostor_metric_level_to = b->get_int("impostor_metric_level_to",impostor_metric_level_to);
        impostor_texture_size = b->get_int("impostor_texture_size",impostor_texture_size);
        leaf_size_mult = b->get_double("leaf_size_mult",leaf_size_mult);
        wood_size_mult = b->get_double("wood_size_mult",wood_size_mult);
        size_diff_factor = b->get_double("size_diff_factor",size_diff_factor);
        size_diff_tolerance = b->get_double("size_diff_tolerance",size_diff_tolerance);
        leaves_opacity = b->get_double("leaves_opacity",leaves_opacity);
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