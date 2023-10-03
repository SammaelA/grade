#pragma once
#include "clustering.h"
#include "helpers.h"
#include "tree_utils/impostor.h"
#include "impostor_metric.h"

class GPUImpostorClusteringHelper : public ImpostorClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
private:
    GLuint results_buf, slices_info_buf, branches_sizes_buf;
};