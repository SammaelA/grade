#include "kmeans_custom.h"
#include "dist_data_table.h"

IntermediateClusteringDataDDT *ddt = nullptr;
//std::pair<float,float> dist_s(int i, int j,)
std::vector<int> KmeansClusteringBase::get_centers_kmeans_plus_plus(
                 IntermediateClusteringData *data, int count)
{
    std::vector<int> centers = {};
    std::vector<BranchClusteringDataImpostor *> possible_centers;
    for (auto *bd : data->branches)
    {
        auto *id = dynamic_cast<BranchClusteringDataImpostor *>(bd);
        if (id && id->can_be_center)
            possible_centers.push_back(id);
    }
    int sz = possible_centers.size();
    count = MIN(count,sz);

    int first_center = urandi(0, sz);
    centers.push_back(first_center);
    std::vector<bool> already_center = std::vector<bool>(sz, false);
    std::vector<double> sums = std::vector<double>(sz+1, 0);
    already_center[first_center] = true;
    for (int i=1;i<count;i++)
    {
        for (int j=0;j<sz;j++)
        {
            double x_sqr = 0;
            if (!already_center[j])
            {
                x_sqr = 1e9;
                for (int &c : centers)
                {
                   // DistData d;
                   // Answer a = ImpostorClusteringHelper::dist_impostor(*(possible_centers[j]),*(possible_centers[c]),
                   //                                                    data->ctx, 0, 1000, &d);
                    auto dist = ddt->ddt.get(j,c);
                    if (dist.first.from < x_sqr)
                        x_sqr = dist.first.from;

                }
            }
            sums[j+1] = sums[j] + SQR(x_sqr);
        }

        double rnd = urand(0, sums[sz-1]);
        for (int j=1;j<sz;j++)
        {
            if (rnd < sums[j])
            {
                centers.push_back(j);
                already_center[j] = true;
                break;
            }
        }
    }
    for (auto &c : centers)
    {
        debug(" %d", c);
    }
    debugnl();
    return centers;
}

bool KmeansClusteringBase::clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result)
{
    ddt = dynamic_cast<IntermediateClusteringDataDDT *>(data);
    data->ctx->self_impostors_raw_atlas = new TextureAtlasRawData(data->ctx->self_impostors_data->atlas);
    
    std::vector<BranchClusteringDataImpostor *> branches;
    for (auto *bd : data->branches)
    {
        auto *id = dynamic_cast<BranchClusteringDataImpostor *>(bd);
        if (id && id->can_be_center)
            branches.push_back(id);
    }
    int sz = branches.size();
    int av_cl_size = get_default_block().get_double("average_cluster_size_goal",10);
    av_cl_size = settings.get_double("average_cluster_size_goal",av_cl_size);
    int clusters_cnt = sz / av_cl_size;
    std::vector<int> centers = get_centers_kmeans_plus_plus(data, clusters_cnt);
    clusters_cnt = centers.size();
    bool change = true;
    int iters = 0;
    std::vector<int> cluster_n = std::vector<int>(sz,-1);
    std::vector<Transform> transfroms = std::vector<Transform>(sz,Transform());
    std::vector<std::pair<int, float>> better_centers = std::vector<std::pair<int, float>>(clusters_cnt, std::pair<int, float>());
    while (change && iters < sz)
    {
        change = false;
        for (int i=0;i<sz;i++)
        {
            float min_d = 1e9;
            float min_r = 0;
            int min_c = -1;
            for (int j=0;j<clusters_cnt;j++)
            {
                int c = centers[j];
                    //DistData d;
                    //Answer a = ImpostorClusteringHelper::dist_impostor(*(branches[i]),*(branches[c]),
                    //                                                   data->ctx, 0, 1000, &d);
                auto dist = ddt->ddt.get(i,c);
                if (dist.first.from < min_d)
                {
                    min_d = dist.first.from;
                    min_r = dist.second.rotation;
                    min_c = j;
                } 
            }
            cluster_n[i] = min_c;
            transfroms[i].rot = min_r;
        }
        
        for (int j=0;j<clusters_cnt;j++)
        {
            better_centers[j].second = 1e20;
        }
        
        for (int i=0;i<sz;i++)
        {
            if (branches[i]->can_be_center)
            {
                float sum = 0;
                for (int j=0;j<sz;j++)
                {
                    if (cluster_n[j] == cluster_n[i])
                    {
                        //DistData d;
                        //Answer a = ImpostorClusteringHelper::dist_impostor(*(branches[i]),*(branches[j]),
                        //                                                   data->ctx, 0, 1000, &d); 
                        auto dist = ddt->ddt.get(i,j);
                        sum += dist.first.from;
                    }
                }
                if (sum < better_centers[cluster_n[i]].second)
                {
                    better_centers[cluster_n[i]].second = sum;
                    better_centers[cluster_n[i]].first = i;
                }
            }
        }
        for (int j=0;j<clusters_cnt;j++)
        {
            if (centers[j] != better_centers[j].first)
            {
                centers[j] = better_centers[j].first;
                change = true;
            }
        }
        iters++;
    }
    delete data->ctx->self_impostors_raw_atlas;

    result = std::vector<ClusterStruct>(clusters_cnt,ClusterStruct());
    for (int j=0;j<clusters_cnt;j++)
    {
        result[j].center = centers[j];
    }
    for (int i=0;i<sz;i++)
    {
        result[cluster_n[i]].members.push_back({i, transfroms[i]});
    }

    return true;
}