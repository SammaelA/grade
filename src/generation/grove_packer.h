#pragma once

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"
#include <unordered_set>
#include "graphics_utils/volumetric_occlusion.h"
#include "tinyEngine/model.h"
#include "graphics_utils/billboard_cloud.h"
#include "core/tree.h"
#include "core/grove.h"
#include "generation/grove_generation_utils.h"
#include "clustering/clustering.h"
#include "save_utils/blk.h"
#include "clustering/clustering.h"
#include "clustering/impostor_metric.h"
#include "generation/generation_settings.h"

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
                            bool visualize_clusters = false, bool save_cluster_data = false);
    static void remove_trees_from_grove(GrovePacked &grove, std::vector<int> &ids);
    void init(Block &packing_params_block);
    void prepare_grove_atlas(GrovePacked &grove, int tex_x, int tex_y, bool save_atlases, bool save_png, 
                             bool alpha_tex_needed);
    GrovePacker() = default;
    explicit GrovePacker(bool shared_ctx);
    ClusteringStrategy get_clustering_strategy() { return cStrategy; }
    std::vector<FullClusteringData *> saved_clustering_data;
    static bool is_valid_tree(::Tree &t);
protected:
    static void recreate_compressed_trees(GrovePacked &grove);
    void add_trees_to_grove_internal(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                                     bool visualize_clusters, bool save_cluster_data);
    void pack_layer(Block &settings, GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                    std::vector<ClusterPackingLayer> &packingLayers, LightVoxelsCube *post_voxels,
                    int layer_from, int layer_to, bool models, bool bill, bool imp,
                    bool visualize_clusters);

    void recalculate_nodes(ClusterData &cl);
    void transform_by_nodes(ClusterData &cl);
    void init();
    void base_init();
    std::vector<ClusterPackingLayer> packingLayersBranches = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrunks = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrees = {ClusterPackingLayer()};

    bool inited = false;
    ClusteringContext *ctx = nullptr;
    ClusteringContext self_ctx;
    Block dummy_block;
    Block settings_block;
    Block *trunks_params = &dummy_block;
    Block *branches_params = &dummy_block;
    Block *trees_params = &dummy_block;
    GroveGenerationData groveGenerationData;
    BranchHeap originalBranches;
    LeafHeap originalLeaves;
    ClusteringStrategy cStrategy = ClusteringStrategy::Merge;
    bool save_clusterizer = false;
    bool shared_context = false;
    int clustering_base_level = 0;

    struct Node
    {
        glm::vec3 position;
    };

    std::map<int, std::vector<Node> > trees_nodes;//for each tree in grove it represents joints of current instance
                                                  //for this tree
};