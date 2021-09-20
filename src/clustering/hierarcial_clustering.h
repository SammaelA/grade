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
    void make_base_clusters(std::vector<BranchClusteringData *> &branches)
    {
        size = branches.size();
        clusters.reserve(size * 2);
        for (int i = 0; i < size; i++)
        {
            clusters.push_back(Cluster(i, branches[i]->can_be_center));
            current_clusters.push_back(i);
            all_clusters_can_be_center = all_clusters_can_be_center && branches[i]->can_be_center;
        }
    }
    void make(int n = 20, int clusters_num = 1);
    Dist get_P_delta(int n, std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta);
    bool all_clusters_can_be_center = true;
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