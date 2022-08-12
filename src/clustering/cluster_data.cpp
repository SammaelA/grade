#include "clustering.h"

bool ClusterData::is_valid()
{
    return (base && !ACDA.clustering_data.empty() &&
            ACDA.clustering_data.size() == IDA.centers_par.size() &&
            ACDA.clustering_data.size() == IDA.centers_self.size() &&
            ACDA.clustering_data.size() == IDA.transforms.size() &&
            ACDA.clustering_data.size() == IDA.type_ids.size() &&
            ACDA.clustering_data.size() == IDA.tree_ids.size());
}