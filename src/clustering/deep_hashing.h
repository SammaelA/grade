#pragma once

#include "../common_utils/python_interaction.h"
#include "hasing.h"

class DeepHashBasedClusteringHelper : public HashBasedClusteringHelper
{
protected:
    BranchClusteringData *convert_branch_deep_hash(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                      BaseBranchClusteringData &data); 
    virtual void branch_conversion_flush(Block &settings, ClusteringContext *ctx) override;

private:
    std::vector<BranchHash *> res_branches;
    std::vector<Branch *> src_branches;
    std::vector<BaseBranchClusteringData> base_branches;
};

class DDTDeepHashBasedClusteringHelper : public DeepHashBasedClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override
    {
        return prepare_intermediate_data_ddt(settings, branches, ctx);
    }
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override
    {
        return convert_branch_deep_hash(settings, base, ctx, data);
    }
};