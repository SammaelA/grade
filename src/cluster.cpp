#include "branch_clusterization.h"

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
            //fprintf(stderr,"ward dist %f %f\n",min,max);
            DistData addData;
            Answer a = Clusterizer::dist(*branch, *(B->branch), min, max, &addData);
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
        //fprintf(stderr,"dist calc %f*%f %f*%f %f*%f \n",d1,a,d2,b,d3,c);
        return distance;
    }
}