#pragma once
#include "common_utils/hash.h"
#include "clustering.h"
#include "hasing.h"
#include "impostor_metric.h"

class KmeansClusteringBase : public ClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override;

private:
    std::vector<int> get_centers_kmeans_plus_plus(IntermediateClusteringData *data, int count);
};