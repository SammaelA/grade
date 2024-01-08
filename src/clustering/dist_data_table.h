#pragma once
#include "clustering.h"
#include "helpers.h"

struct DistDataTable;

struct DistData
{
    float dist = 1;
    float rotation = 0;
    DistData(float _dist = 1, float _rotation = 0)
    {
        dist = _dist;
        rotation = _rotation;
    }
};
struct Cluster
{
    int branch_n = -1;
    Cluster *U = nullptr;
    Cluster *V = nullptr;
    int size = 0;
    bool can_be_center = true;
    std::map<Cluster *, DistData> distances;
    Cluster(int n, bool can_be_center);
    Cluster(Cluster *_U, Cluster *_V);
    void to_base_clusters(std::vector<Cluster *> &clusters);
    float ward_dist(Cluster *B, float min = 1.0, float max = 0.0);
};
struct DistDataTable
{
    std::pair<Answer, DistData> *data = nullptr;
    int n;
    void create(int _n)
    {
        n = _n;
        data = safe_new<std::pair<Answer, DistData>>(n * n, "DistDataTable_data");
    }
    void clear()
    {
        safe_delete<std::pair<Answer, DistData>>(data, "DistDataTable_data");
    }
    inline std::pair<Answer, DistData> get(int x, int y) { return data[x * n + y]; }
    inline void set(int x, int y, Answer &a, DistData &d)
    {
        data[x * n + y] = std::pair<Answer, DistData>(a, d);
    }
    int size() { return n; }
    void copy(DistDataTable &ddt)
    {
        n = ddt.n;
        data = safe_new<std::pair<Answer, DistData>>(n * n, "DistDataTable_data");
        for (int i=0;i<n*n;i++)
          data[i] = ddt.data[i];
    }
    ~DistDataTable() { clear(); }
};

struct IntermediateClusteringDataDDT : public IntermediateClusteringData
{
    DistDataTable ddt;

    virtual void clear(){};
};