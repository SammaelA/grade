#pragma once
#include "clustering.h"
#include "impostor_metric.h"

struct ClusterPackingLayer;
void visualize_clusters(Block &settings, std::vector<BranchClusteringData *> branches, 
                        std::vector<ClusteringBase::ClusterStruct> &clusters,
                        ClusteringContext *ctx, std::string file_name, int w, int h);

void prepare_dataset(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers);
void save_csv(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers);