#pragma once

#include "hasing.h"
class ImpostorHashClusteringHelper2 : public HashBasedClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
    //BranchClusteringData that is returned after convert_branch could be not fully initialized.
    //call this function before using it.  
    virtual void branch_conversion_flush(Block &settings, ClusteringContext *ctx) override;
private:
    std::vector<BranchHash *> res_branches;
    std::vector<Branch *> src_branches;
    std::vector<BaseBranchClusteringData> base_branches;
};