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
#include "generation/generation_task.h"

class Heightmap;
struct ClusterAdditionalData
{
    std::list<InstancedBranch>::iterator instanced_branch;
    std::vector<std::list<BillboardData>::iterator> small_billboards;
    std::vector<std::list<BillboardData>::iterator> large_billboards;
    std::list<Impostor>::iterator impostor;
    bool is_presented = false;
    bool has_instanced_branch = false;
    bool has_impostor = false;
};
struct ClusterPackingLayer
{
    std::vector<ClusterData> clusters;
    std::vector<ClusterAdditionalData> additional_data;
};
struct GrovePackingParams
{
  GrovePackingParams(unsigned _task, const ImpostorBaker::ImpostorGenerationParams &igp)
  {
    task = _task;
    impostor_generation_params = igp;
  }
  unsigned task = MINIMUM_FOR_RENDER;
  ImpostorBaker::ImpostorGenerationParams impostor_generation_params;
};
class GrovePacker
{
public:
    void add_trees_to_grove(const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                            bool visualize_clusters = false, bool save_cluster_data = false);
    void remove_trees(GrovePacked &grove, std::vector<int> &ids);
    void init(const Block &packing_params_block, const std::vector<TreeTypeData> &types);
    void prepare_grove_atlas(GrovePacked &grove, int tex_x, int tex_y, bool save_atlases, bool save_png, 
                             bool alpha_tex_needed);
    GrovePacker() = default;
    explicit GrovePacker(bool shared_ctx);
    ClusteringStrategy get_clustering_strategy() { return cStrategy; }
    std::vector<FullClusteringData *> saved_clustering_data;
    static bool is_valid_tree(::Tree &t);
protected:
    static void recreate_compressed_trees(GrovePacked &grove);
    void add_trees_to_grove_internal(const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                                     bool visualize_clusters, bool save_cluster_data);
    void pack_layer(Block &settings, const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                    std::vector<ClusterPackingLayer> &packingLayers,
                    int layer_from, int layer_to, bool models, bool bill, bool imp,
                    bool visualize_clusters);

    void recalculate_nodes(ClusterData &cl);
    void transform_by_nodes(ClusterData &cl);
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
    std::vector<TreeTypeData> types;
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