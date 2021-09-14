#pragma once
#include "clustering.h"
#include "hasing.h"

enum XKmeans
{
    K_MEANS,
    X_MEANS,
};
class XKmeansPyClusteringBase : public ClusteringBase
{
protected:
    virtual bool clusterize_internal(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result, XKmeans xkmeans,
                                     float initial_clusters_frac, bool force_recalculate_centers);
};
class KmeansPyClusteringBase : public XKmeansPyClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override
    {
        return clusterize_internal(settings, data, result, XKmeans::K_MEANS, 1, false);
    }
};
class XmeansPyClusteringBase : public XKmeansPyClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override
    {
        return clusterize_internal(settings, data, result, XKmeans::X_MEANS, 0.5, true);
    }
};


enum Optics_DBscan
{
    OPTICS,
    DBSCAN,
};
class OpticsDBscanPyClusteringBase : public ClusteringBase
{
protected:
    virtual bool clusterize_internal(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result, 
                                     Optics_DBscan type);
};
class DBscanPyClusteringBase : public OpticsDBscanPyClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override
    {
        return clusterize_internal(settings, data, result, Optics_DBscan::DBSCAN);
    }
};
class OpticsPyClusteringBase : public OpticsDBscanPyClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override
    {
        return clusterize_internal(settings, data, result, Optics_DBscan::OPTICS);
    }
};