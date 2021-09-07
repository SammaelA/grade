#pragma once
#include "clustering.h"
#include "helpers.h"
#include "../impostor.h"

struct BranchClusteringDataImpostor : public BranchClusteringData
{
    std::list<Impostor>::iterator self_impostor;
    virtual void clear() override;
};

/*struct ImpostorClusteringContext : public ClusteringContext
{
    ImpostorsData *self_impostors_data = nullptr;
    TextureAtlasRawData *self_impostors_raw_atlas = nullptr;
    virtual void clear() override;
    ~ImpostorClusteringContext() { clear(); }
};*/

class ImpostorClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override;
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
};