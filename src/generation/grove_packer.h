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
#include "common_utils/blk.h"
#include "clustering/clustering.h"
#include "clustering/impostor_metric.h"
#include "generation/generation_task.h"

class Heightmap;
struct ClusterAdditionalData
{
    friend class boost::serialization::access;

    std::list<InstancedBranch>::iterator instanced_branch;//points in grove.instancedBranches list

    std::vector<std::list<BillboardData>::iterator> small_billboards;//unused right now
    std::vector<std::list<BillboardData>::iterator> large_billboards;//unused right now
    std::list<Impostor>::iterator impostor;//unused right now
    bool is_presented = false;
    bool has_instanced_branch = false;
    bool has_impostor = false;

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if (!Archive::is_loading::value)
      {
        int pos = cur_ser_helper->instanced_branches.pos_by_it(instanced_branch);
        ar & pos;
      }
      else
      {
        int pos = -1;
        ar & pos;
        instanced_branch = cur_ser_helper->instanced_branches.it_by_pos(pos);
      }
      ar & is_presented;
      ar & has_instanced_branch;
      ar & has_impostor;
    }
};
struct ClusterPackingLayer
{
    friend class boost::serialization::access;

    std::vector<ClusterData> clusters;
    std::vector<ClusterAdditionalData> additional_data;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & clusters;
      ar & additional_data;
    }
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
    friend class boost::serialization::access;

    void add_trees_to_grove(const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                            bool visualize_clusters = false, bool save_cluster_data = false);
    void remove_trees(GrovePacked &grove, std::vector<int> &ids);
    void init(const Block &packing_params_block, const std::vector<TreeTypeData> &types);
    void prepare_grove_atlas(GrovePacked &grove, int tex_x, int tex_y, bool save_atlases, bool save_png, 
                             bool alpha_tex_needed);
    void post_serialization_recover(GrovePacked &grove);
    void clear();
    ~GrovePacker() {clear();}
    GrovePacker() = default;
    GrovePacker& operator=(GrovePacker&&) = delete;
    ClusteringStrategy get_clustering_strategy() { return cStrategy; }
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
    std::list<InstancedBranch>::iterator pack_cluster(ClusterData &cluster, GrovePacked &grove, int lvl_from, int lvl_to);

    std::vector<ClusterPackingLayer> packingLayersBranches = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrunks = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrees = {ClusterPackingLayer()};
    bool inited = false;
    ClusteringContext ctx;
    Block settings_block;
    Block trunks_params;
    Block branches_params;
    Block trees_params;
    BranchHeap originalBranches;
    LeafHeap originalLeaves;
    ClusteringStrategy cStrategy = ClusteringStrategy::Merge;
    bool save_clusterizer = false;
    int clustering_base_level = 0;
    std::map<int, std::vector<glm::vec3> > trees_nodes;//for each tree in grove it represents joints of current instance
                                                       //for this tree
    int ib_id_counter = 1;
private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & inited;
      ar & ctx;
      ar & packingLayersBranches;
      ar & packingLayersTrees;
      ar & packingLayersTrunks;
      ar & settings_block;
      ar & branches_params;
      ar & trees_params;
      ar & trunks_params;
      ar & originalLeaves;
      ar & originalBranches;
      ar & cStrategy;
      ar & save_clusterizer;
      ar & clustering_base_level;
      ar & trees_nodes;

      //we need to save and restore this global id counters to prevent id collision
      ar & ib_id_counter;
      ar & cur_cluster_id;
    }
};