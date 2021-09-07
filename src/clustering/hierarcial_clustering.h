#pragma once
#include "clustering.h"
#include "helpers.h"
#include "dist_data_table.h"

struct ClusterDendrogramm
{
    struct Dist
    {
        int U;
        int V;
        float d;
        Dist(int _U, int _V, float _d)
        {
            U = _U;
            V = _V;
            d = _d;
        }
    };
    int size;
    Cluster *root;
    std::vector<Cluster> clusters;
    std::list<int> current_clusters;
    void make_base_clusters(int count)
    {
        size = count;
        clusters.reserve(count * 2);
        for (int i = 0; i < count; i++)
        {
            clusters.push_back(Cluster(i));
            current_clusters.push_back(i);
        }
    }
    void make(int n = 20, int clusters_num = 1);
    Dist get_P_delta(int n, std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta);
};

class HierarcialClusteringBase : public ClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override;
    friend class Cluster;

private:
    int get_typical(std::vector<Cluster *> &clusters);
    static Answer get_dist(int n1, int n2, DistData *data);
    IntermediateClusteringDataDDT *main_data;
};