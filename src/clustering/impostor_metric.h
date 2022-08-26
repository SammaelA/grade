#pragma once
#include "clustering.h"
#include "helpers.h"
#include "graphics_utils/impostor.h"
#include "impostor_similarity_params.h"

struct BranchClusteringDataImpostor : public BranchClusteringData
{
    friend class boost::serialization::access;

    std::list<Impostor>::iterator self_impostor;
    virtual void clear() override;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<BranchClusteringData>(*this);
      if (!Archive::is_loading::value)
      {
        int pos = cur_ser_helper->clust_metric_impostors.pos_by_it(self_impostor);
        ar & pos;
      }
      else
      {
        int pos = -1;
        ar & pos;
        self_impostor = cur_ser_helper->clust_metric_impostors.it_by_pos(pos);
      }
    }
};

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