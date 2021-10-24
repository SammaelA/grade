#pragma once
#include "clustering.h"
#include "helpers.h"
#include "../impostor.h"
#include "impostor_similarity_params.h"

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
struct DistData;
class ImpostorClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override;
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override;
    static Answer dist_impostor(BranchClusteringDataImpostor &bwd1, BranchClusteringDataImpostor &bwd2, 
                                ClusteringContext *current_data, float min, float max, DistData *data);
};

extern ImpostorSimilarityParams isimParams;