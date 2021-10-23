#include "clustering_benchmark.h"
#include <chrono>
#include "dist_data_table.h"
#include "clustering_debug_utils.h"

#define STEPS (ClusteringStep::CLUSTERING_STEPS_COUNT)
struct ClusteringResult
{
    std::string name;
    float total_time_ms;
    float average_cluster_size[STEPS];
    float average_in_cluster_distance[STEPS];//the distance is calculated based on REFERENCE DDT
};
struct ReferenceResult
{
    
};
void print_clustering_result(ClusteringResult &res, bool is_reference, ClusteringResult &reference)
{
    debug("--------------------\n");
    if (is_reference)
    {
        debug("Name: %s (reference)\n",res.name.c_str());
        debug("Time spent %f seconds\n",1e-4*res.total_time_ms);
        debug("Average cluster size %.1f/%.1f/%.1f\n",
              res.average_cluster_size[0], res.average_cluster_size[1], res.average_cluster_size[2]);
        debug("Average clustering error size %.3f/%.3f/%.3f\n",
              res.average_in_cluster_distance[0], res.average_in_cluster_distance[1], res.average_in_cluster_distance[2]);
    }
    else
    {
        debug("Name: %s\n",res.name.c_str());
        debug("Time spent %f seconds %d %% from reference\n",1e-4*res.total_time_ms, 
              (int)(100*res.total_time_ms/reference.total_time_ms));
        debug("Average cluster size %.1f/%.1f/%.1f  %.1f%%/%.1f%%/%.1f%% from reference\n",
              res.average_cluster_size[0], res.average_cluster_size[1], res.average_cluster_size[2],
              (100*res.average_cluster_size[0]/reference.average_cluster_size[0]),
              (100*res.average_cluster_size[1]/reference.average_cluster_size[1]),
              (100*res.average_cluster_size[2]/reference.average_cluster_size[2]));
        debug("Average clustering error size %.3f/%.3f/%.3f  %.1f%%/%.1f%%/%.1f%% from reference\n",
              res.average_in_cluster_distance[0], res.average_in_cluster_distance[1], res.average_in_cluster_distance[2],
              (100*res.average_in_cluster_distance[0]/reference.average_in_cluster_distance[0]),
              (100*res.average_in_cluster_distance[1]/reference.average_in_cluster_distance[1]),
              (100*res.average_in_cluster_distance[2]/reference.average_in_cluster_distance[2]));
    }
    debug("--------------------\n");
}
void ClusteringBenchmark::perform_benchmark(std::string benchmark_blk_path, AbstractTreeGenerator *gen, 
                                            GroveGenerationData &ggd, Heightmap *h)
{
    debug("starting clustering benchmark. Preparing grove.\n");
    bool need_visualize_clusters = true;
    int d = debug_level;
    debug_level = (debug_level == 11) ? 11 : 1000;
    float generation_time = 0;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    ::Tree *trees = new ::Tree[ggd.trees_count];
    gen->create_grove(ggd, trees, *h);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    generation_time  = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("Original grove created. Took %f seconds\n",1e-4*generation_time);
    
    BlkManager man;
    Block settings;
    Block default_settings;
    man.load_block_from_file(benchmark_blk_path, settings);
    std::vector<ClusteringResult> results;
    settings.add_block("default", &default_settings);


    std::map<int,int> pos_in_table_by_id[STEPS];
    DistDataTable DDTs[STEPS];
    FullClusteringData *clusts[STEPS];

    for (int i=settings.size() - 1;i>=0;i--)
    {
        bool reference = (i == settings.size() - 1);
        Block *bl = settings.get_block(i);
        if (!bl)
        {
            logerr("clustering benchmark settings block should contain only blocks with clustering settings");
            continue;
        }
        results.emplace_back();
        results.back().name = settings.get_name(i);
        debug("starting clustering with settings %s\n",results.back().name.c_str());
        GrovePackerStat packer(true);
        groves.emplace_back();
        
        t1 = std::chrono::steady_clock::now();
        packer.init(*bl);
        
        if (reference)
            packer.start_save_clusterizer();
        packer.add_trees_to_grove(ggd, groves.back(),trees,h, false, true);
        if (reference && packer.saved_clustering_data.size() == STEPS)
        {
            for (int j=0;j<STEPS;j++)
            {
                auto *clust = packer.saved_clustering_data[j];
                if (!clust)
                {
                    logerr("empty clust");
                    continue;
                }
                IntermediateClusteringDataDDT *ddt_data = dynamic_cast<IntermediateClusteringDataDDT *>(clust->id);
                if (!ddt_data)
                {
                    logerr("reference clustring algorithm don't have DDT in intermediate data");
                }
                DDTs[j] = ddt_data->ddt;
                pos_in_table_by_id[j] = clust->pos_in_table_by_branch_id;
                clusts[j] = clust;
                for (auto &pair : pos_in_table_by_id[j])
                {
                    debugl(11,"pair %d-->%d\n",pair.first,pair.second);
                }
                for (int x = 0; x < MIN(10,DDTs[j].n);x++)
                {
                    for (int y = 0; y < MIN(10,DDTs[j].n);y++)
                    {
                        debugl(11,"%.1f ",DDTs[j].get(x,y).first.from);
                    }
                    debugl(11,"\n");
                }
            }
        }
        t2 = std::chrono::steady_clock::now();
        results.back().total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        
        for (int j=0;j<STEPS;j++)
        {
            auto structure = packer.clusterStructures((ClusteringStep)j);
            int total_branches = 0;
            float total_dist = 0;
            for (auto &str : structure)
            {
                debugl(11,"cluster %d{",str.base_id);
                total_branches += str.ids.size();
                int base_pos = pos_in_table_by_id[j].at(str.base_id);
                for (int id : str.ids)
                {
                    debugl(11,"%d,",id);
                    int pos = pos_in_table_by_id[j].at(id);
                    total_dist += CLAMP(DDTs[j].get(base_pos,pos).first.from,0,1);
                    debugl(11,"(%d %d %f)",base_pos,pos,CLAMP(DDTs[j].get(base_pos,pos).first.from,0,1));
                }
                debugl(11,"}\n");
            }
            debugl(11,"\n");
            results.back().average_in_cluster_distance[j] = total_dist/total_branches;
            results.back().average_cluster_size[j] = (float)total_branches/(structure.size());
            
            if (need_visualize_clusters && ((ClusteringStep)j == ClusteringStep::BRANCHES))
            {
              auto *clust = packer.saved_clustering_data[j];
              if (!clust)
              {
                logerr("empty clust");
                continue;
              }
              std::string vis_name = results.back().name + "_clusters_" + std::to_string(j);
              auto icon_sizes = clusts[j]->ctx->self_impostors_data->atlas.get_slice_size();
              visualize_clusters(settings, clusts[j]->id->branches, clust->clusters, clusts[j]->ctx, vis_name,icon_sizes.x,icon_sizes.y);
            }
        }
    }

    debug("Clustering benchmark finished %d algoritms compared\n",results.size());

    for (int i=0;i<results.size();i++)
    {
        print_clustering_result(results[i],i==0,results[0]);
    }
    debug_level = d;
}

std::vector<ClusterPackingLayer> &GrovePackerStat::clusterLayers(ClusteringStep step)
{
    if (step == ClusteringStep::TRUNKS)
        return packingLayersTrunks;
    else if (step == ClusteringStep::BRANCHES)
        return packingLayersBranches;
    else if (step == ClusteringStep::TREES)
        return packingLayersTrees;
}
std::vector<GrovePackerStat::ClusterStructure> GrovePackerStat::clusterStructures(ClusteringStep step)
{
    std::vector<ClusterPackingLayer> &layers = clusterLayers(step);
    std::vector<GrovePackerStat::ClusterStructure> structures;
    for (ClusterPackingLayer &layer : layers)
    {
        for (auto &cd : layer.clusters)
        {
            structures.emplace_back();
            structures.back().base_id = cd.base->self_id;
            structures.back().ids = cd.ACDA.ids;
        }
    }

    return structures;
}