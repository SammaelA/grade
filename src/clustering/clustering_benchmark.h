#pragma once
#include "../grove_packer.h"
#include "../terrain.h"
#include "../abstract_generator.h"

class ClusteringBenchmark
{
public:
    void perform_benchmark(std::string benchmark_blk_path, AbstractTreeGenerator *gen, GroveGenerationData &ggd, Heightmap *h);
    GrovePacked &get_grove(int n) {return groves[CLAMP(n,0,groves.size()-1)];}
    int grove_count() {return groves.size();}
private:
    std::vector<GrovePacked> groves;
};

class GrovePackerStat : public GrovePacker
{
public:
    struct ClusterStructure
    {
        int base_id;
        std::vector<int> ids;
    };
    void start_save_clusterizer() { save_clusterizer = true;}
    std::vector<ClusterPackingLayer> &clusterLayers(ClusteringStep step);
    std::vector<ClusterStructure> clusterStructures(ClusteringStep step);
};