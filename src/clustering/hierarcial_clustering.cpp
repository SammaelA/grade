#include "hierarcial_clustering.h"

struct HierarcialClusteringParams
{
    float average_cluster_size_goal = 4;

    void load_from_block(Block *b);
    void load(Block *b);
};
void HierarcialClusteringParams::load(Block *b)
{
    if (!b)
        return;
    
    Block &def = get_default_block();
    load_from_block(&def);
    load_from_block(b);
}
void HierarcialClusteringParams::load_from_block(Block *b)
{
    average_cluster_size_goal = b->get_double("average_cluster_size_goal",average_cluster_size_goal);
}

HierarcialClusteringParams hsParams;
DistDataTable *currentDdt = nullptr;

bool HierarcialClusteringBase::clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result)
{
    main_data = dynamic_cast<IntermediateClusteringDataDDT *>(data);
    hsParams.load(&settings);
    if (!main_data)
    {
        logerr("Error: HierarcialClusteringBase got intermediate data of a wrong class");
        return false;
    }
    if (main_data->elements_count <= 0)
    {
        logerr("Error: HierarcialClusteringBase got empty intermediate data.");
        return false;
    }
    if (main_data->elements_count != main_data->branches.size() || main_data->elements_count != main_data->ddt.size())
    {
        logerr("Error: HierarcialClusteringBase size mismatch in intermediate data");
        return false;
    }

    currentDdt = &(main_data->ddt);
    ClusterDendrogramm den;
    den.make_base_clusters(main_data->branches);
    den.make();

    for (int c_num : den.current_clusters)
    {
        Cluster &c = den.clusters[c_num];
        std::vector<Cluster *> base_clusters;
        c.to_base_clusters(base_clusters);

        result.emplace_back();
        int center = get_typical(base_clusters);
        result.back().center = center;
        for (auto *cl : base_clusters)
        {
            auto p = currentDdt->get(center,cl->branch_n);
            Transform t;
            t.rot = p.second.rotation;
            result.back().members.push_back({cl->branch_n,t});
        }
    }
}

ClusterDendrogramm::Dist
ClusterDendrogramm::get_P_delta(int n, std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta)
{
    debugl(1, "get P delta\n");
    Dist md(-1, -1, 1e9);
    int k = n > current_clusters.size() ? current_clusters.size() : n;
    k = n > current_clusters.size() / 3 ? n : current_clusters.size() / 3;
    int n1 = 0;
    auto i = current_clusters.begin();
    auto j = current_clusters.rbegin();
    while (n1 < k)
    {
        if (*i == *j)
        {
        }
        else
        {
            float distance = clusters[*i].ward_dist(&(clusters[*j]), md.d);
            if (distance < md.d && distance > 0.001)
            {
                md.d = distance;
                md.U = *i;
                md.V = *j;
            }
        }
        i++;
        j++;
        n1++;
    }
    delta = md.d;
    debugl(1, "P delta delta %f\n", delta);
    for (int u : current_clusters)
    {
        for (int v : current_clusters)
        {
            if (u != v)
            {
                float distance = clusters[u].ward_dist(&(clusters[v]), delta);
                if (distance <= delta)
                {
                    P_delta.push_back(Dist(u, v, distance));
                    if (distance < md.d)
                    {
                        md.d = distance;
                        md.U = u;
                        md.V = v;
                    }
            }
            }
        }
    }
    debugl(1,"P_delta size = %d md = (%d %d %f)\n", P_delta.size(), md.U, md.V, md.d);
    return md;
}
void ClusterDendrogramm::make(int n, int clusters_num)
{
    int initial_clusters = current_clusters.size();

    if (!all_clusters_can_be_center)
    {
        //first we need to make clusters with at least one branch that
        //can be center in each
        std::list<int> no_centers_clusters;
        std::list<int> centers_clusters;
        for (int cl : current_clusters)
        {
            if (!clusters[cl].can_be_center)
                no_centers_clusters.push_back(cl);
            else
                centers_clusters.push_back(cl);
        }

        for (int cl : no_centers_clusters)
        {
            float min = 1000;
            int min_pos = -1;
            for (int cent : centers_clusters)
            {
                float d = clusters[cl].ward_dist(&(clusters[cent]));
                if (d < min)
                {
                    min = d;
                    min_pos = cent;
                }
            }
            current_clusters.remove(cl);
            current_clusters.remove(min_pos);
            clusters.push_back(Cluster(&(clusters[cl]),&(clusters[min_pos])));
            int W = clusters.size() - 1;
            current_clusters.push_back(W);

            centers_clusters.remove(min_pos);
            centers_clusters.push_back(W);
        }
    }

    std::list<Dist> P_delta;
    float delta;
    Dist min = get_P_delta(n, current_clusters, P_delta, delta);
    for (int i = 1; i < size; i++)
    {
        if (P_delta.empty())
            min = get_P_delta(n, current_clusters, P_delta, delta);
        else if (min.U < 0)
        {
            for (Dist &d : P_delta)
            {
                if (d.U == d.V)
                    debugl(0, "error in P_delta %d\n", d.U);
                if (d.d < min.d)
                {
                    min.U = d.U;
                    min.V = d.V;
                    min.d = d.d;
                }
            }
        }
        if (min.d >= 1000 || current_clusters.size() <= clusters_num)
        {
            break;
            //makes no sense to merge clusters with maximum distance between them.
        }
        else if ((float)initial_clusters/current_clusters.size() >= hsParams.average_cluster_size_goal &&
                 hsParams.average_cluster_size_goal > 0)
        {
            break;
        }
        current_clusters.remove(min.U);
        current_clusters.remove(min.V);
        clusters.push_back(Cluster(&(clusters[min.U]), &(clusters[min.V])));
        int W = clusters.size() - 1;
        auto dit = P_delta.begin();
        while (dit != P_delta.end())
        {
            Dist &d = *dit;
            if (d.U == min.U || d.V == min.U || d.U == min.V || d.V == min.V)
                dit = P_delta.erase(dit);
            else
                dit++;
        }
        for (int S : current_clusters)
        {
            float d = clusters[W].ward_dist(&(clusters[S]), delta);
            if (d < delta)
            {
                P_delta.push_back(Dist(W, S, d));
            }
        }
        current_clusters.push_back(W);
        min = Dist(-1, -1, 1000);
    }
    int sum = 0;
    for (int S : current_clusters)
    {
        debugl(1, "cluster %d size = %d\n", S, clusters[S].size);
        sum += clusters[S].size;
    }
    float comp = (float)size/current_clusters.size();
    debugl(17, "%d %d %f clusters elements compression\n", current_clusters.size(), size, comp);
}

Answer HierarcialClusteringBase::get_dist(int n1, int n2, DistData *data)
{
    auto p = currentDdt->get(n1,n2);
    if (data)
        *data = p.second;
    return p.first;
}
int HierarcialClusteringBase::get_typical(std::vector<Cluster *> &clusters)
{
    float min_d = 1e9;
    int min_pos = 0;
    float d_0;
    for (int i = 0; i<clusters.size();i++)
    {
        if (!main_data->branches[clusters[i]->branch_n]->can_be_center)
            continue;
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
    if (!main_data->branches[clusters[min_pos]->branch_n]->can_be_center)
    {
        logerr("H clustering error - created a cluster without possible centers!");
    }
    return clusters[min_pos]->branch_n;
}
Cluster::Cluster(int n, bool _can_be_center)
{
    branch_n = n;
    size = 1;
    can_be_center = _can_be_center;
}
Cluster::Cluster(Cluster *_U, Cluster *_V)
{
    U = _U;
    V = _V;
    size = U->size + V->size;
    can_be_center = _U->can_be_center || _V->can_be_center;
}
void Cluster::to_base_clusters(std::vector<Cluster *> &clusters)
{
    if (branch_n >= 0)
        clusters.push_back(this);
    else if (U && V)
    {
        U->to_base_clusters(clusters);
        V->to_base_clusters(clusters);
    }
}
float Cluster::ward_dist(Cluster *B, float min, float max)
{
    if (this == B)
        return 0;
    auto it = distances.find(B);
    if (it != distances.end())
        return it->second.dist;
    it = B->distances.find(this);
    if (it != B->distances.end())
        return it->second.dist;
    if (branch_n >= 0)
    { //dist between single branch and cluster. Reduce second cluster
        if (B->branch_n >= 0)
        { //dist between single branches.
            DistData addData;
            Answer a = HierarcialClusteringBase::get_dist(branch_n, B->branch_n, &addData);
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