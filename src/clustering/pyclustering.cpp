#include "pyclustering.h"
#include "../../dependencies/pyclustering/ccore/include/pyclustering/definitions.hpp"
#include "../../dependencies/pyclustering/ccore/include/pyclustering/cluster/kmeans.hpp"
#include "../../dependencies/pyclustering/ccore/include/pyclustering/cluster/xmeans.hpp"
#include "../../dependencies/pyclustering/ccore/include/pyclustering/cluster/kmeans_plus_plus.hpp"
double L2(const std::vector<double> &v1,const std::vector<double> &v2)
{
    double res = 0;
    for (int i=0;i<MIN(v1.size(),v2.size());i++)
    {
        res += SQR(v1[i] - v2[i]);
    }

    return res;
}

bool XKmeansPyClusteringBase::clusterize_internal(Block &settings, IntermediateClusteringData *data, 
                                                  std::vector<ClusterStruct> &result,
                                                  XKmeans xkmeans, float initial_clusters_frac, bool force_recalculate_centers)
{
    IntermediateClusteringDataVectorsList *main_data = dynamic_cast<IntermediateClusteringDataVectorsList *>(data);

    if (!main_data)
    {
        logerr("Error: KmeansPyClusteringBase got intermediate data of a wrong class");
        return false;
    }
    if (main_data->branches.size() != main_data->feature_vectors.size())
    {
        logerr("Error: KmeansPyClusteringBase got incorrect intermediate data");
        return false;
    }

    float average_cluster_size_goal = get_default_block().get_double("average_cluster_size_goal",10);
    average_cluster_size_goal = settings.get_double("average_cluster_size_goal",average_cluster_size_goal);
    bool recalculate_centers = get_default_block().get_bool("kmeans_recalculate_centers",true);
    recalculate_centers = settings.get_bool("kmeans_recalculate_centers",recalculate_centers);
    recalculate_centers = recalculate_centers && force_recalculate_centers;

    int max_centers = main_data->feature_vectors.size()/average_cluster_size_goal;
    int clusters = MAX(max_centers*initial_clusters_frac,1);
    
    pyclustering::dataset centers;
    pyclustering::clst::index_sequence centers_ns;
    if (data->all_branches_can_be_centers)
    {
        pyclustering::clst::kmeans_plus_plus centers_search(clusters,main_data->feature_vectors.size());
        centers_search.initialize(main_data->feature_vectors,centers_ns);
    }
    else
    {
        pyclustering::dataset possible_centers;
        std::map<int, int> index_by_pos_center_index;
        for (int i=0;i<main_data->feature_vectors.size();i++)
        {
            if (main_data->branches[i]->can_be_center)
            {
                possible_centers.push_back(main_data->feature_vectors[i]);
                index_by_pos_center_index.emplace(possible_centers.size() - 1, i);
            }
        }
        clusters = MIN(clusters, possible_centers.size());

        pyclustering::clst::kmeans_plus_plus centers_search(clusters,possible_centers.size());
        pyclustering::clst::index_sequence tmp_centers_ns;
        centers_search.initialize(possible_centers,centers_ns);
        for (auto id : tmp_centers_ns)
        {
            centers_ns.emplace_back(index_by_pos_center_index.at(id));
        }
    }
    for (auto id : centers_ns)
    {
        centers.push_back(main_data->feature_vectors[id]);
    }
    pyclustering::clst::kmeans_data kmeans_data;
    if (xkmeans == K_MEANS)
    {
        pyclustering::clst::kmeans kmeans(centers);
        kmeans.process(main_data->feature_vectors, kmeans_data);
    }
    else if (xkmeans == X_MEANS)
    {
        pyclustering::clst::xmeans xmeans(centers, max_centers - clusters);
        xmeans.process(main_data->feature_vectors, kmeans_data);
    }

    if (centers_ns.size() != kmeans_data.clusters().size() && !force_recalculate_centers)
    {
        logerr("centers count changed and all centers should be recalculated");
        recalculate_centers = true;
    }
    for (int i=0;i<kmeans_data.clusters().size();i++)
    {
        const auto &cluster = kmeans_data.clusters()[i];
        const auto &center = kmeans_data.centers()[i];
        result.emplace_back();
        if (!recalculate_centers)
        {
            for (auto &id : cluster)
            {
                ClusteringBase::Transform t;
                t.rot = 0;
                result.back().members.push_back({(int)id, t});
            }
            result.back().center = centers_ns[i];
        }
        else
        {
            float min_dist = 1e9;
            int best_center = -1;
            for (auto &id : cluster)
            {
                ClusteringBase::Transform t;
                t.rot = 0;
                result.back().members.push_back({(int)id, t});
                if (main_data->branches[id]->can_be_center)
                {
                    float dist = L2(center, main_data->feature_vectors[id]);
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        best_center = id;
                    }
                }
            }
            if (best_center == -1)
            {
                logerr("clustering finished with cluster where no member can be a center. TODO: deal with it");
            }
            else
            {
                result.back().center = best_center;
            }
        }
    }

    return true;
}