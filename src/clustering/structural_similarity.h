#pragma once
#include "branch_with_data.h"
#include "dist_data_table.h"

class StructuralSimilarityClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override;
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override = 0;
};
class GPUSSClusteringHelper : public StructuralSimilarityClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
};
class CPUSSClusteringHelper : public StructuralSimilarityClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
};