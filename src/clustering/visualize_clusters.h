#pragma once
#include "clustering.h"
#include "impostor_metric.h"

void visualize_clusters(Block &settings, std::vector<BranchClusteringData *> branches, 
                        std::vector<ClusteringBase::ClusterStruct> &clusters,
                        ClusteringContext *ctx, std::string file_name, int w, int h);
