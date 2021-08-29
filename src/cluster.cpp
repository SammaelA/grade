#include "branch_clusterization.h"
Clusterizer *Clusterizer::Cluster::currentClusterizer = nullptr;
Clusterizer::Cluster::Cluster(BranchWithData *bwd)
{
    branch = bwd;
    size = 1;
}
Clusterizer::Cluster::Cluster(Cluster *_U, Cluster *_V)
{
    U = _U;
    V = _V;
    size = U->size + V->size;
}
void Clusterizer::Cluster::to_branch_data(std::vector<Branch *> &branches)
{
    if (branch)
        branches.push_back(branch->b);
    else if (U && V)
    {
        U->to_branch_data(branches);
        V->to_branch_data(branches);
    }
}
void Clusterizer::Cluster::to_base_clusters(std::vector<Cluster *> &clusters)
{
    if (branch)
        clusters.push_back(this);
    else if (U && V)
    {
        U->to_base_clusters(clusters);
        V->to_base_clusters(clusters);
    }
}
float Clusterizer::Cluster::ward_dist(Cluster *B, float min, float max)
{
    auto it = distances.find(B);
    if (it != distances.end())
        return it->second.dist;
    it = B->distances.find(this);
    if (it != B->distances.end())
        return it->second.dist;
    if (branch)
    { //dist between single branch and cluster. Reduce second cluster
        if (B->branch)
        { //dist between single branches.
            DistData addData;
            Answer a = currentClusterizer->get_dist(*branch, *(B->branch), &addData);
            float distance = a.to;
            addData.dist = distance;
            if (a.exact)
                distances.emplace(B, addData);
            return distance;
        }
        else
        {
            float d1 = B->U->ward_dist(this);
            float d2 = B->V->ward_dist(this);
            float d3 = B->U->ward_dist(B->V);
            float a = (float)(size + B->U->size) / (size + B->size);
            float b = (float)(size + B->V->size) / (size + B->size);
            float c = -(float)size / (size + B->size);
            float distance = a * d1 + b * d2 + c * d3;
            distances.emplace(B, distance);
            return distance;
        }
    }
    else
    { //dist between two clusters. Reduce first cluster
        float d1 = U->ward_dist(B);
        float d2 = V->ward_dist(B);
        float d3 = U->ward_dist(V);
        float a = (float)(B->size + U->size) / (size + B->size);
        float b = (float)(B->size + V->size) / (size + B->size);
        float c = -(float)(B->size) / (size + B->size);
        float distance = a * d1 + b * d2 + c * d3;
        distances.emplace(B, distance);
        return distance;
    }
}
bool valid(ClusterData &cd)
{
    return (cd.base && !cd.ACDA.originals.empty() &&
            cd.ACDA.originals.size() == cd.ACDA.rotations.size() &&
            cd.ACDA.originals.size() == cd.IDA.centers_par.size() &&
            cd.ACDA.originals.size() == cd.IDA.centers_self.size() &&
            cd.ACDA.originals.size() == cd.IDA.transforms.size() &&
            cd.ACDA.originals.size() == cd.IDA.type_ids.size() &&
            cd.ACDA.originals.size() == cd.IDA.tree_ids.size());
}
Branch *Clusterizer::Cluster::prepare_to_replace(std::vector<ClusterData> &base_clusters, InstanceDataArrays &IDA, 
                                                 AdditionalClusterDataArrays &ADCA, long &cluster_id)
{
    std::vector<Cluster *> clusters;
    return prepare_to_replace(base_clusters, IDA, ADCA, clusters, cluster_id);
}
Branch *Clusterizer::Cluster::prepare_to_replace(std::vector<ClusterData> &base_clusters, InstanceDataArrays &IDA, 
                                                 AdditionalClusterDataArrays &ADCA, std::vector<Cluster *> &clusters,
                                                 long &cluster_id)
{
    to_base_clusters(clusters);
    if (clusters.empty())
        return nullptr;
    BranchWithData *br = get_typical(clusters);
    glm::mat4 base_transform_inv = glm::inverse(br->transform);
    for (Cluster *cl : clusters)
    {
        glm::mat4 tr = (cl->branch->transform) * base_transform_inv;
        if (cl->branch->base_cluster_id < 0 || cl->branch->base_cluster_id >= base_clusters.size())
        {
            logerr("failed to merge clusters. Base cluster is not found in given array and will be ignored");
        }
        else
        {
            ClusterData &base = base_clusters[cl->branch->base_cluster_id];
            if (!valid(base))
            {
                logerr("failed to merge clusters. Base cluster is not valid and will be ignored");
            }
            else
            {
                int sz = base.ACDA.originals.size();
                for (int i=0;i<sz;i++)
                {
                    IDA.transforms.push_back(base.IDA.transforms[i]*tr);
                    IDA.centers_par.push_back(base.IDA.centers_par[i]);
                    IDA.centers_self.push_back(base.IDA.centers_self[i]);
                    IDA.type_ids.push_back(base.IDA.type_ids[i]);
                    IDA.tree_ids.push_back(base.IDA.tree_ids[i]);
                    ADCA.rotations.push_back(base.ACDA.rotations[i]);
                    ADCA.originals.push_back(base.ACDA.originals[i]);
                    ADCA.ids.push_back(base.ACDA.ids[i]);
                }
            }
            if (br == cl->branch)
                cluster_id = base.id;
        }
    }
    return br->original;
}
Clusterizer::BranchWithData *Clusterizer::Cluster::get_typical(std::vector<Cluster *> &clusters)
{
    float min_d = 1e9;
    int min_pos = 0;
    float d_0;
    for (int i = 0; i<clusters.size();i++)
    {
        float d = 0;
        for (int j = 0; j<clusters.size();j++)
        {
            if (i==j)
                continue;
            d += clusters[i]->ward_dist(clusters[j]);
        }
        if (d < min_d)
        {
            min_d = d;
            min_pos = i;
        }
        if (i == 0)
            d_0 = d;
    }

    return (clusters[min_pos])->branch;
}