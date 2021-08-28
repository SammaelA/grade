#include "clustering_benchmark.h"
#include <chrono>

struct ClusteringResult
{
    std::string name;
    float total_time_ms;
};
void print_clustering_result(ClusteringResult &res, bool is_reference, ClusteringResult &reference)
{
    debug("--------------------\n");
    if (is_reference)
    {
        debug("Name: %s (reference)\n",res.name.c_str());
        debug("Time spent %f seconds\n",1e-4*res.total_time_ms);
    }
    else
    {
        debug("Name: %s\n",res.name.c_str());
        debug("Time spent %f seconds %d %% from reference\n",1e-4*res.total_time_ms, 
              (int)(100*res.total_time_ms/reference.total_time_ms));
    }
    debug("--------------------\n");
}
void ClusteringBenchmark::perform_benchmark(std::string benchmark_blk_path, AbstractTreeGenerator *gen, 
                                            GroveGenerationData &ggd, Heightmap *h)
{
    debug("starting clustering benchmark. Preparing grove.\n");
    int d = debug_level;
    debug_level = 1000;
    float generation_time = 0;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    ::Tree *trees = new ::Tree[ggd.trees_count];
    gen->create_grove(ggd, trees, *h);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    generation_time  = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    BlkManager man;
    Block settings;
    man.load_block_from_file(benchmark_blk_path, settings);
    std::vector<ClusteringResult> results;
    debug("Original grove created. Took %f seconds\n",1e-4*generation_time);
    for (int i=0;i<settings.size();i++)
    {
        Block *bl = settings.get_block(i);
        if (!bl)
        {
            logerr("clustering benchmark settings block should contain only blocks with clustering settings");
            continue;
        }
        results.emplace_back();
        results.back().name = settings.get_name(i);
        debug("starting clustering with settings %s\n",results.back().name.c_str());
        GrovePacker packer;
        groves.emplace_back();
        t1 = std::chrono::steady_clock::now();
        packer.init(*bl);
        packer.add_trees_to_grove(ggd, groves.back(),trees,h);
        t2 = std::chrono::steady_clock::now();
        results.back().total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    }

    debug("Clustering benchmark finished %d algoritms compared\n",results.size());

    for (int i=0;i<results.size();i++)
    {
        print_clustering_result(results[i],i==0,results[0]);
    }
    debug_level = d;
}