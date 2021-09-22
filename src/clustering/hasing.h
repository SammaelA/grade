#pragma once
#include "clustering.h"
#include "helpers.h"
#include "../hash.h"

struct BranchHash : public BranchClusteringData
{
    std::vector<Hash> hashes;
};

class HashBasedClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) = 0;
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) = 0;
protected:
    BranchClusteringData *convert_branch_impostor(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                  BaseBranchClusteringData &data);
    BranchClusteringData *convert_branch_eigin_vectors(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                       BaseBranchClusteringData &data);
    BranchClusteringData *convert_branch_impostor_dct(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                      BaseBranchClusteringData &data);
    IntermediateClusteringData *prepare_intermediate_data_ddt(Block &settings, std::vector<BranchClusteringData *> branches,
                                                              ClusteringContext *ctx);
    IntermediateClusteringData *prepare_intermediate_data_simple(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                 ClusteringContext *ctx);
};

class DDTHashBasedClusteringHelper : public HashBasedClusteringHelper
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
        return convert_branch_eigin_vectors(settings, base, ctx, data);
    }
};


struct IntermediateClusteringDataVectorsList : public IntermediateClusteringData
{
    std::vector<std::vector<double>> feature_vectors;
};

class SimpleHashBasedClusteringHelper : public HashBasedClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override
    {
        return prepare_intermediate_data_simple(settings, branches, ctx);
    }
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override
    {
        return convert_branch_eigin_vectors(settings, base, ctx, data);
    }
};

class DDTImpostorHashClusteringHelper : public HashBasedClusteringHelper
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
        return convert_branch_impostor(settings, base, ctx, data);
    }
};

class SimpleImpostorHashClusteringHelper : public HashBasedClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override
    {
        return prepare_intermediate_data_simple(settings, branches, ctx);
    }
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override
    {
        return convert_branch_impostor(settings, base, ctx, data);
    }
};

class DDTImpostorDCTHashClusteringHelper : public HashBasedClusteringHelper
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
        return convert_branch_impostor_dct(settings, base, ctx, data);
    }
};

class SimpleImpostorDCTHashClusteringHelper : public HashBasedClusteringHelper
{
public:
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override
    {
        return prepare_intermediate_data_simple(settings, branches, ctx);
    }
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override
    {
        return convert_branch_impostor_dct(settings, base, ctx, data);
    }
};