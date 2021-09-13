#pragma once
#include "clustering.h"
#include "helpers.h"
#include "../impostor.h"
#include "impostor_metric.h"

class GPUImpostorClusteringHelper : public ImpostorClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
private:
    struct Task
    {
        glm::uvec4 from;
        glm::uvec4 to;
    };
    GLuint tasks_buf, results_buf, raw_data_buf;
};