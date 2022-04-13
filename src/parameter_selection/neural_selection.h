#pragma once
#include "core/tree.h"
#include "core/scene.h"
#include "common_utils/parameter.h"
#include "graphics_utils/texture_atlas.h"

class NeuralEstimator
{
public:
    void prepare_dataset(ParameterList &params, GroveGenerationData &tree_ggd, 
                         std::string save_path, int impostor_size = 256, int impostor_slices = 8,
                         int batch_size = 1024, int batches_cnt = 256, bool save_metainfo = true,
                         int joints_limit = 100000);
private:
};