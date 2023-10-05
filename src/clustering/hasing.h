#pragma once
#include "clustering.h"
#include "helpers.h"
#include "common_utils/hash.h"

struct BranchHash : public BranchClusteringData
{
  friend class boost::serialization::access;

  std::vector<Hash> hashes;

private:
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    logerr("BranchHash cannot be serialized. Serialization is not implemented");
  }
};
struct ImpostorSimilarityParams;
class HashBasedClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                 BaseBranchClusteringData &data) override = 0;
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override;
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override = 0;
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
    void create_impostor_temp(Block &settings, Branch *base, ClusteringContext *ctx, BaseBranchClusteringData &data,
                              ImpostorSimilarityParams &isimParams, TextureAtlas &temp_atlas, Impostor &imp);
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