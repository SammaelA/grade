#pragma once

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"
#include <unordered_set>
#include "volumetric_occlusion.h"
#include "tinyEngine/utility/model.h"
#include "billboard_cloud.h"
#include "tree.h"
#include "grove.h"
#include "grove_generation_utils.h"
#include "clustering/clustering.h"
#include "tinyEngine/save_utils/blk.h"
#include "clustering/clustering.h"
#include "clustering/impostor_metric.h"

class DebugVisualizer;
class Heightmap;
struct ClusterAdditionalData
{
    std::list<InstancedBranch>::iterator instanced_branch;
    std::vector<std::list<BillboardData>::iterator> small_billboards;
    std::vector<std::list<BillboardData>::iterator> large_billboards;
    std::list<Impostor>::iterator impostors;
    bool is_presented = false;
};
struct ClusterPackingLayer
{
    std::vector<ClusterData> clusters;
    std::vector<ClusterAdditionalData> additional_data;
};

class GrovePacker
{
public:
    void add_trees_to_grove(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                            bool visualize_clusters = false);
    void init(Block &packing_params_block);
    std::vector<FullClusteringData *> saved_clustering_data;
protected:
    void add_trees_to_grove_internal(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                                     bool visualize_clusters);
    void pack_layer(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                std::vector<ClusterPackingLayer> &packingLayers, LightVoxelsCube *post_voxels,
                int layer_from, int layer_to, bool models, bool bill, bool imp,
                bool visualize_clusters);
    void transform_all_according_to_root(GrovePacked &grove);
    void init();
    std::vector<ClusterPackingLayer> packingLayersBranches = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrunks = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrees = {ClusterPackingLayer()};

    bool inited = false;
    ClusteringContext ctx;
    Block settings_block;
    GroveGenerationData groveGenerationData;
    BranchHeap originalBranches;
    LeafHeap originalLeaves;

    bool save_clusterizer = false;
};