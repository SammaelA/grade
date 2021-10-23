#pragma once

#include <string>

struct ClusteringDebugInfo
{
    bool prepare_dataset = false;
    std::string dataset_name = "dataset";

    bool save_csv = false;
    std::string csv_file_name = "clusters.csv";

    bool visualize_clusters = false;
    std::string visualize_clusters_name = "clusters";

    bool progress_bars = false;
};

extern ClusteringDebugInfo clusteringDebugInfo;