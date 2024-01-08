#include "common_utils/utility.h"
#include "clustering/clustering.h"
#include "tinyEngine/engine.h"
#include "tree_utils/tree_modeling.h"
#include "common_utils/distribution.h"
#include <math.h>
#include <algorithm>
#include "common_utils/body.h"
#include <chrono>
#include "tree_utils/impostor.h"
#include "graphics_utils/terrain.h"
#include "common_utils/field_2d.h"
#include "generation/grove_generation_utils.h"
#include "synthetic_trees_generator.h"
#include "generation/grove_packer.h"
#include "clustering/default_clustering_params.h"
#include "clustering/clustering_debug_utils.h"
#include "clustering/clustering_debug_status.h"
#include "metainfo_manager.h"
#include "scene_generator_helper.h"
#include <set>

void pack_branch_recursively(::Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, int lvl_from, int lvl_to);
std::list<InstancedBranch>::iterator GrovePacker::pack_cluster(ClusterData &cluster, GrovePacked &grove, 
                                                               int lvl_from, int lvl_to)
{
    grove.instancedBranches.push_back(InstancedBranch());
    grove.instancedBranches.back().IDA = cluster.IDA;
    grove.instancedBranches.back().id = ib_id_counter++;
    std::vector<unsigned> &ids = grove.instancedBranches.back().branches;
    ::Branch *base = cluster.base;
    pack_branch_recursively(base, grove, ids, lvl_from, lvl_to);
    grove.instancedBranches.back().bbox = BillboardCloudRaw::get_minimal_bbox(base);

    std::list<InstancedBranch>::iterator it = grove.instancedBranches.end();
    it--;
    return it;
}

void pack_branch_recursively(::Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, int lvl_from, int lvl_to)
{
    if (b->level > lvl_to)
        return;
    if (b->level >= lvl_from)
    {
        PackedBranch pb;
        b->pack(pb);
        ids.push_back(grove.instancedCatalogue.add(pb, b->level));
    }
    for (::Joint &j : b->joints)
    {
        for (::Branch *br : j.childBranches)
            pack_branch_recursively(br, grove, ids, lvl_from, lvl_to);
    }
}
void add_occluder(LightVoxelsCube *voxels, Branch *b)
{
    for (Joint &j : b->joints)
    {
        voxels->set_occluder(j.pos, j.leaf ? 2.5 : 1);
        for (Branch *chb : j.childBranches)
            add_occluder(voxels, chb);
    }
}
void add_occluder(LightVoxelsCube *voxels, Tree *trees, int count)
{
    for (int i = 0; i < count; i++)
    {
        if (trees[i].valid && trees[i].root)
            add_occluder(voxels, trees[i].root);
    }
}
bool GrovePacker::is_valid_tree(::Tree &t)
{
    if (!t.valid)
    {
        //logerr("tree %u was marked as not valid by a generator. Grove packer will ignore it", t.id);
        return false;
    }
    if (!(t.valid && t.leaves && t.branchHeaps.size() >= 2))
    {
        //logerr("tree %u do not have some of essential data structures. %d Grove packer will ignore it", t.id,
        //       t.branchHeaps.size());
        return false;
    }
    if (t.branchHeaps[0]->branches.empty())
    {
        //logerr("tree %u do not have trunk. Grove packer will ignore it", t.id);
        return false;
    }
    if (t.branchHeaps[1]->branches.empty())
    {
        //logerr("tree %u do not have branches instead of trunk. Grove packer will ignore it", t.id);
        return false;
    }
    return true;
}

Clusterizer2 *get_clusterizer(Block &settings, ClusteringStrategy &cStrategy)
{
  std::string c_strategy_name = "merge";
  c_strategy_name = settings.get_string("clustering_strategy", c_strategy_name);

  if (c_strategy_name == "recreate")
  {
    cStrategy = ClusteringStrategy::Recreate;
  }
  else
  {
    cStrategy = ClusteringStrategy::Merge;
  }
  Clusterizer2 *cl = new Clusterizer2(cStrategy);
  cl->prepare(settings);

  return cl;
}

void GrovePacker::pack_layer(Block &settings, const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                             std::vector<ClusterPackingLayer> &packingLayers,
                             int layer_from, int layer_to, bool models, bool bill, bool imp,
                             bool visualize_clusters)
{
    unsigned task = params.task;

    bool need_something = task & (GenerationTask::SYNTS | GenerationTask::CLUSTERIZE) ||
                          ((task & GenerationTask::MODELS) && models) ||
                          ((task & GenerationTask::BILLBOARDS) && bill) ||
                          ((task & GenerationTask::IMPOSTORS) && imp);
    if (!need_something)
      return;

    int clusters_before = packingLayers[clustering_base_level].clusters.size();
    Clusterizer2 *cl = get_clusterizer(settings, cStrategy);

    struct ClusterInfo
    {
        int layer;
        int pos;
    };
    struct ExpandInfo
    {
        ClusterInfo from;
        ClusterInfo to;
    };
    std::vector<ClusterInfo> new_clusters;
    std::vector<ClusterInfo> removed_clusters;
    std::vector<ExpandInfo> expanded_clusters;
    std::vector<int> cleared_layers;

    /* && cStrategy == ClusteringStrategy::Merge*/
    if ((task & (GenerationTask::CLUSTERIZE)))
    {
        cl->get_base_clusters(settings, trees_external, trees_count, layer_from, packingLayers[clustering_base_level].clusters, &ctx, true);
    }
    else
    {
        //just add base clusters to zero level cluster list
        cl->get_base_clusters(settings, trees_external, trees_count, layer_from, packingLayers[clustering_base_level].clusters, 
                              &ctx, false);
    }

    //we do not modify previously existed clusters on this step. Only add some new
    for (int i = clusters_before; i < packingLayers[clustering_base_level].clusters.size(); i++)
    {
        packingLayers[clustering_base_level].additional_data.emplace_back();
        new_clusters.push_back(ClusterInfo{clustering_base_level, i});
    }

    if (task & (GenerationTask::CLUSTERIZE))
    {
        int max_clusters_in_layer = settings.get_int("max_clusters_in_layer",1000);
        int max_layer = settings.get_int("recursive_clustering_layers",1);
        int max_lookup_level = (cStrategy == ClusteringStrategy::Merge ? packingLayers.size() - 1 : clustering_base_level);
        max_lookup_level = packingLayers.size() - 1;
        for (int i = clustering_base_level; i <= max_lookup_level; i++)
        {
            int prev_size = packingLayers[i].clusters.size();
            if (i < max_layer && prev_size > max_clusters_in_layer)
            {
                cleared_layers.push_back(i);
                if (i + 1 == packingLayers.size())
                {
                    packingLayers.push_back(ClusterPackingLayer());
                }
                std::map<long, int> old_cl_poses;
                BitVector remains;
                remains.resize(prev_size, false);
                int old_size = packingLayers[i + 1].clusters.size();
                cl->clusterize(settings, packingLayers[i].clusters, packingLayers[i + 1].clusters, &ctx,
                               save_clusterizer, visualize_clusters);
                int new_size = packingLayers[i + 1].clusters.size();

                for (int j = 0; j < prev_size; j++)
                {
                    old_cl_poses.emplace(packingLayers[i].clusters[j].id, j);
                    debugl(6, "old cluster %d %d\n", j, (int)(packingLayers[i].clusters[j].id));
                }
                for (int j = old_size; j < new_size; j++)
                {
                    auto it = old_cl_poses.find(packingLayers[i + 1].clusters[j].id);
                    if (it == old_cl_poses.end())
                    {
                        //it is a new cluster
                        new_clusters.push_back(ClusterInfo{i + 1, j});
                    }
                    else
                    {
                        //this cluster is an extension of existing one
                        remains.set(it->second, true);
                        expanded_clusters.push_back(ExpandInfo{ClusterInfo{i, it->second}, ClusterInfo{i + 1, j}});
                    }
                    debugl(6, "new cluster %d %d\n", j, (int)(packingLayers[i + 1].clusters[j].id));
                    packingLayers[i + 1].additional_data.emplace_back();
                }
                for (int j = 0; j < prev_size; j++)
                {
                    if (!remains.get(j))
                    {
                        //old cluster should be deleted
                        removed_clusters.push_back(ClusterInfo{i, j});
                    }
                }
            }
        }
    }
    bool first = true;

    for (auto &info : new_clusters)
    {
        Branch *original = packingLayers[info.layer].clusters[info.pos].base;
        Branch *new_original = originalBranches.new_branch();

        for (auto &info : expanded_clusters)
        {
            if (packingLayers[info.from.layer].clusters[info.from.pos].base == original)
            {
                packingLayers[info.from.layer].clusters[info.from.pos].base = new_original;
            }
            if (packingLayers[info.to.layer].clusters[info.to.pos].base == original)
            {
                packingLayers[info.to.layer].clusters[info.to.pos].base = new_original;
            }
        }
        new_original->deep_copy(original,originalBranches, &originalLeaves);
        packingLayers[info.layer].clusters[info.pos].base = new_original;
      
        if (layer_from == 0 && layer_to == 0)
            recalculate_nodes(packingLayers[info.layer].clusters[info.pos]); 
        else if (layer_from == 1)
            transform_by_nodes(packingLayers[info.layer].clusters[info.pos]);
    
        if (models && (task & (GenerationTask::MODELS)))
        {
            auto it = pack_cluster(packingLayers[info.layer].clusters[info.pos], grove, layer_from, layer_to);
            packingLayers[info.layer].additional_data[info.pos].instanced_branch = it;
            packingLayers[info.layer].additional_data[info.pos].has_instanced_branch = true;
            packingLayers[info.layer].additional_data[info.pos].is_presented = true;
        }
        if (bill && (task & (GenerationTask::BILLBOARDS)))
        {
            BillboardCloudRaw cloud;
            cloud.prepare(ULTRALOW, layer_from, packingLayers[info.layer].clusters[info.pos], ctx.types,
                          &grove.clouds[2], packingLayers[info.layer].additional_data[info.pos].large_billboards);
            cloud.prepare(ULTRALOW, layer_from + 1, packingLayers[info.layer].clusters[info.pos], ctx.types,
                          &grove.clouds[3], packingLayers[info.layer].additional_data[info.pos].small_billboards);
        }
        if (imp && (task & (GenerationTask::IMPOSTOR_FULL_GROVE)))
        {
            debugl(6, "full grove impostors are temporary unavailable\n");
        }
        if (imp && (task & (GenerationTask::IMPOSTORS)))
        {
          packingLayers[info.layer].additional_data[info.pos].has_impostor = true;
            ImpostorBaker ib;
            ib.prepare(params.impostor_generation_params, layer_from, packingLayers[info.layer].clusters[info.pos], 
                       ctx.types, &(grove.impostors[1]), packingLayers[info.layer].additional_data[info.pos].impostor,
                    new_clusters.size());
        }

        debugl(6, "added cluster %d\n", packingLayers[info.layer].clusters[info.pos].id);
    }
    for (auto &info : expanded_clusters)
    {       
        if (layer_from == 0 && layer_to == 0)
            recalculate_nodes(packingLayers[info.to.layer].clusters[info.to.pos]); 
        else if (layer_from == 1)
            transform_by_nodes(packingLayers[info.to.layer].clusters[info.to.pos]);
        
        if (models && (task & (GenerationTask::MODELS)))
        {
            auto &it = packingLayers[info.from.layer].additional_data[info.from.pos].instanced_branch;
            it->IDA = packingLayers[info.to.layer].clusters[info.to.pos].IDA;
            packingLayers[info.to.layer].additional_data[info.to.pos].instanced_branch = it;
            packingLayers[info.to.layer].additional_data[info.to.pos].has_instanced_branch = true;
            packingLayers[info.to.layer].additional_data[info.to.pos].is_presented = true;
        }
        if (bill && (task & (GenerationTask::BILLBOARDS)))
        {
            BillboardCloudRaw cloud;
            auto &from_a = packingLayers[info.from.layer].additional_data[info.from.pos];
            auto &to_a = packingLayers[info.to.layer].additional_data[info.to.pos];
            to_a.large_billboards = from_a.large_billboards;
            to_a.small_billboards = from_a.small_billboards;

            cloud.extend(ULTRALOW, layer_from, packingLayers[info.to.layer].clusters[info.to.pos], ctx.types,
                         &grove.clouds[2], to_a.large_billboards);
            cloud.extend(ULTRALOW, layer_from + 1, packingLayers[info.to.layer].clusters[info.to.pos], ctx.types,
                         &grove.clouds[3], to_a.small_billboards);

            packingLayers[info.to.layer].additional_data[info.to.pos].is_presented = true;
        }
        if (imp && (task & (GenerationTask::IMPOSTORS)))
        {
            auto &from_a = packingLayers[info.from.layer].additional_data[info.from.pos];
            auto &to_a = packingLayers[info.to.layer].additional_data[info.to.pos];

            to_a.has_impostor = from_a.has_impostor;
            to_a.impostor = from_a.impostor;
            to_a.impostor->IDA = packingLayers[info.to.layer].clusters[info.to.pos].IDA;
        }
        debugl(6, "replaced cluster %d with %d\n", packingLayers[info.from.layer].clusters[info.from.pos].id,
               packingLayers[info.to.layer].clusters[info.to.pos].id);
    }
    for (auto &info : removed_clusters)
    {
        bool need_clear_data = true;
        if (need_clear_data)
        {
            if (models && (task & (GenerationTask::MODELS)))
            {
                auto &it = packingLayers[info.layer].additional_data[info.pos].instanced_branch;
                for (unsigned &id : it->branches)
                    grove.instancedCatalogue.remove(id);
                grove.instancedBranches.erase(it);
            }
            if (bill && (task & (GenerationTask::BILLBOARDS)))
            {
                for (auto &it : packingLayers[info.layer].additional_data[info.pos].large_billboards)
                {
                    for (auto &bill : it->billboards)
                    {
                        grove.clouds[2].atlas.remove_tex(bill.id);
                    }
                    grove.clouds[2].billboards.erase(it);
                }

                for (auto &it : packingLayers[info.layer].additional_data[info.pos].small_billboards)
                {
                    for (auto &bill : it->billboards)
                    {
                        grove.clouds[3].atlas.remove_tex(bill.id);
                    }
                    grove.clouds[3].billboards.erase(it);
                }
            }
            if (imp && (task & (GenerationTask::IMPOSTORS)))
            {
                auto &it = packingLayers[info.layer].additional_data[info.pos].impostor;
                for (auto &bill : it->slices)
                {
                    grove.impostors[1].atlas.remove_tex(bill.id);
                }
                grove.impostors[1].atlas.remove_tex(it->top_slice.id);
                grove.impostors[1].impostors.erase(it);
            }

            packingLayers[info.layer].clusters[info.pos].base->mark_dead();
            debugl(6, "removed cluster %d\n", packingLayers[info.layer].clusters[info.pos].id);
        }
    }
    for (auto i : cleared_layers)
    {
        packingLayers[i].additional_data.clear();
        packingLayers[i].clusters.clear();
        debugl(6, "level cleared %d\n", i);
    }
    int total_cl_cnt = 0;
  int br_cl_cnt = 0;
  std::vector<ClusterPackingLayer *> cpls;
  for (auto &cpl : packingLayersBranches)
    cpls.push_back(&cpl);
  for (auto &cpl : packingLayersTrunks)
    cpls.push_back(&cpl);
  for (auto &cpl : packingLayersTrees)
    cpls.push_back(&cpl);
  for (int i=0;i<cpls.size();i++)
  {
    if (cpls[i]->clusters.size() != cpls[i]->additional_data.size())
          logerr("Clusted packing layer is corrupted. %d clusters, %d additional data", cpls[i]->clusters.size(), cpls[i]->additional_data.size());
    for (int j=0;j<cpls[i]->clusters.size();j++)
    {
      total_cl_cnt++;
      if (cpls[i]->additional_data[j].has_instanced_branch)
      {
        br_cl_cnt++;
      }
    }
  }
  if (br_cl_cnt != grove.instancedBranches.size())
    logerr("we have %d clusters, %d with branches and %d instanced branches",total_cl_cnt, br_cl_cnt, grove.instancedBranches.size());

  if (cl)
    delete cl;
}

void GrovePacker::recalculate_nodes(ClusterData &cl)
{
    for (int k = 0; k < cl.IDA.transforms.size(); k++)
    {
        auto nit = trees_nodes.find(cl.IDA.tree_ids[k]);
        if (nit == trees_nodes.end())
        {
            logerr("cannot find list of nodes for tree %d", cl.IDA.tree_ids[k]);
        }
        else
        {
            nit->second.clear();
            for (Joint &j : cl.base->joints)
            {
                glm::vec3 pos = glm::vec3(cl.IDA.transforms[k] * glm::vec4(j.pos, 1));
                nit->second.push_back({pos});
            }
        }
    }
}

void GrovePacker::transform_by_nodes(ClusterData &cl)
{
    glm::vec4 base_pos = glm::vec4(cl.base->joints.front().pos, 1);
    for (int i = 0; i < cl.IDA.transforms.size(); i++)
    {
        glm::vec3 pos = glm::vec3(cl.IDA.transforms[i] * base_pos);
        float min_dist = 1000;
        glm::vec3 min_vec = glm::vec3(0, 0, 0);
        auto it = trees_nodes.find(cl.IDA.tree_ids[i]);
        if (it == trees_nodes.end())
        {
            logerr("cannot find list of nodes for tree %d", cl.IDA.tree_ids[i]);
        }
        else
        {
            for (auto &node : it->second)
            {
                glm::vec3 dist = node - pos;
                float d = glm::dot(dist, dist);
                if (d < min_dist)
                {
                    min_dist = d;
                    min_vec = dist;
                }
            }
            cl.IDA.transforms[i] = glm::translate(glm::mat4(1.0f), min_vec) * cl.IDA.transforms[i];
        }
    }
}

void GrovePacker::add_trees_to_grove(const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                                     bool visualize_clusters, bool save_cluster_data)
{
    int max_trees_in_patch = settings_block.get_int("max_trees_in_patch",1000);
    visualize_clusters = visualize_clusters || clusteringDebugInfo.visualize_clusters;
    if (visualize_clusters || clusteringDebugInfo.prepare_dataset || clusteringDebugInfo.save_csv)
    {
        add_trees_to_grove_internal(params, grove, trees_external, trees_count, visualize_clusters, true);

        if (clusteringDebugInfo.prepare_dataset)
        {
            bool status = prepare_directory(clusteringDebugInfo.dataset_name);
            if (status)
                prepare_dataset(clusteringDebugInfo.dataset_name, &ctx, packingLayersBranches);
            else
                logerr("unable to create directory to save dataset. Exiting.");
        }
        if (clusteringDebugInfo.save_csv)
            save_csv(clusteringDebugInfo.csv_file_name, &ctx, packingLayersBranches);
    }
    else if (trees_count <= max_trees_in_patch)
    {
        add_trees_to_grove_internal(params, grove, trees_external, trees_count, visualize_clusters, save_cluster_data);
    }
    else
    {
        int start_pos = 0;
        int t_cnt = MIN(max_trees_in_patch, trees_count - start_pos);
        while (start_pos < trees_count)
        {
            add_trees_to_grove_internal(params, grove, trees_external + start_pos, t_cnt, visualize_clusters, save_cluster_data);
            start_pos+=max_trees_in_patch;
        }
    }
}

bool clear_and_validate_branch_structure(Branch *br)
{
  if (br->joints.size() < 2)
    return false; 
  else if (br->joints.size() != br->segments.size() + 1)
    return false;
  for (Joint &j : br->joints)
  {
    auto it = j.childBranches.begin();
    while (it != j.childBranches.end())
    {
      bool valid = clear_and_validate_branch_structure(*it);
      if (valid)
        it++;
      else
      {
        (*it)->joints.clear();
        (*it)->segments.clear();
        it = j.childBranches.erase(it);
      }
    }
  }
  return true;
}

bool clear_and_validate_tree_structure(::Tree &tree)
{
  return (tree.root->level == 0) && clear_and_validate_branch_structure(tree.root);
}

void GrovePacker::add_trees_to_grove_internal(const GrovePackingParams &params, GrovePacked &grove, ::Tree *trees_external, int trees_count,
                                              bool visualize_clusters, bool save_cluster_data)
{
    if (!inited)
    {
      logerr("Grove Packer not inited");
      return;
    }
    
    for (int i = grove.impostors.size(); i < 2; i++)
    {
        grove.impostors.push_back(ImpostorsData());
        grove.impostors.back().valid = true;
    }
    for (int i = grove.clouds.size(); i < 4; i++)
        grove.clouds.push_back(BillboardCloudData());

    save_clusterizer = save_cluster_data;
    int valid_trees_cnt = 0;
    Branch dummy_branch;
    for (int i = 0; i < trees_count; i++)
    {
        bool valid = is_valid_tree(trees_external[i]);
        if (valid)
          valid = clear_and_validate_tree_structure(trees_external[i]);
        if (!valid)
        {
            //dummy tree
            trees_external[i].clear();
            trees_external[i] = Tree();
            trees_external[i].root = &dummy_branch;
        }
        trees_external[i].valid = valid;
        valid_trees_cnt += valid;
    }
    if (trees_count == 0)
    {
        logerr("Grove is empty. ");
        return;
    }
    if (valid_trees_cnt == 0)
    {
        logerr("Grove has %d tree(s) but none of them are valid.", trees_count);
        return;
    }

    debugl(1, "Packing grove with %d/%d valid trees\n", valid_trees_cnt, trees_count);

    for (int i = 0; i < trees_count; i++)
    {
        auto it = trees_nodes.emplace(trees_external[i].id,std::vector<glm::vec3>{});
        for (Joint &j : trees_external[i].root->joints)
        {
            it.first->second.push_back(j.pos);
        }
    }


    current_clustering_step = ClusteringStep::TRUNKS;
    pack_layer(trunks_params, params, grove, trees_external, trees_count, packingLayersTrunks,
               0, 0, true, false, false, false);

    current_clustering_step = ClusteringStep::BRANCHES;
    pack_layer(branches_params, params, grove, trees_external, trees_count, packingLayersBranches,
               1, 1000, true, true, false, visualize_clusters);
    if (params.task & (GenerationTask::IMPOSTORS))
    {
      //TODO: support impostors
      current_clustering_step = ClusteringStep::TREES;
      pack_layer(trees_params, params, grove, trees_external, trees_count, packingLayersTrees,
                0, 1000, false, false, true, false);
    }
    recreate_compressed_trees(grove);

    originalBranches.clear_removed();
    originalLeaves.clear_removed();
}

void GrovePacker::recreate_compressed_trees(GrovePacked &grove)
{

  grove.instancedBranchesDirect.clear();
  grove.instancedBranchesDirect.reserve(grove.instancedBranches.size());
  for (auto it = grove.instancedBranches.begin(); it != grove.instancedBranches.end(); it++)
    grove.instancedBranchesDirect.emplace_back(it);

  // fill tree structures
  // TODO: work only with tree structures that changed after clustering
  std::set<int> p_t_ids;
  for (auto &t : grove.compressedTrees)
    p_t_ids.emplace(t.global_id);
  grove.compressedTrees.clear();
  grove.trees_by_global_id.clear();
  int j = 0;
  for (auto &ib : grove.instancedBranches)
  {
    int branch_level = 1000;
    if (ib.branches.empty())
    {
      logerr("instancedBranches contains malformed branch");
      continue;
    }
    else
    {
      branch_level = grove.instancedCatalogue.get(ib.branches[0]).level;
      if (branch_level > 1)
      {
        logerr("instancedBranches with level > 1 are not supported correctly for branch structure creation. It will be added to trunk directly");
      }
    }
    for (int i = 0; i < ib.IDA.centers_par.size(); i++)
    {
      auto it = grove.trees_by_global_id.find(ib.IDA.tree_ids[i]);
      if (it == grove.trees_by_global_id.end())
      {
        int global_id = ib.IDA.tree_ids[i];

        it = grove.trees_by_global_id.emplace(global_id, grove.compressedTrees.size()).first;
        grove.compressedTrees.emplace_back();
        grove.compressedTrees.back().global_id = ib.IDA.tree_ids[i];
        grove.compressedTrees.back().LOD_roots.emplace_back();
      }
      auto &root = grove.compressedTrees[it->second].LOD_roots[0];
      if (branch_level == 0 && ib.branches.size() == 1)//it's trunk
      {
        if (root.model_num >= 0)
          logerr("tree %d has more than one root node. It's strange", ib.IDA.tree_ids[i]);
        else
        {
          grove.compressedTrees[it->second].pos = ib.IDA.centers_self[i];
          root.model_num = j;
          root.instance_num = i;
        }
      }
      else
      {
        root.children.push_back(CompressedTree::Node(CompressedTree::MODEL, j, i));
      }
    }
    j++;
  }

  for (auto &t : grove.compressedTrees)
  {
    glm::vec3 min_pos = glm::vec3(1e9,1e9,1e9);
    glm::vec3 max_pos = glm::vec3(-1e9,-1e9,-1e9);
    std::vector<CompressedTree::Node> nodes = t.LOD_roots;
    while (!nodes.empty())
    {
      std::vector<CompressedTree::Node> new_nodes;
      for (auto &node : nodes)
      {
        if (node.type == CompressedTree::MODEL)
        {
          auto &bb = grove.instancedBranchesDirect[node.model_num]->bbox;
          auto &tr = grove.instancedBranchesDirect[node.model_num]->IDA.transforms[node.instance_num];
          glm::mat4 rot_inv(glm::vec4(bb.a, 0), glm::vec4(bb.b, 0), glm::vec4(bb.c, 0), glm::vec4(0, 0, 0, 1));
          glm::vec3 pos = rot_inv * glm::vec4(bb.position, 1.0f);
          glm::vec3 p1 = tr * glm::vec4(pos,1);
          glm::vec3 p2 = tr * glm::vec4(pos + bb.sizes.x*bb.a + bb.sizes.y*bb.b + bb.sizes.z*bb.c,1);
          min_pos = min(min_pos, p1);
          min_pos = min(min_pos, p2);
          max_pos = max(max_pos, p1);
          max_pos = max(max_pos, p2);
        }
        new_nodes.insert(new_nodes.end(), node.children.begin(), node.children.end());
      }
      nodes = std::move(new_nodes);
    }
    t.bbox = AABB(min_pos, max_pos);
  }
}

void GrovePacker::init(const Block &packing_params_block, const std::vector<TreeTypeData> &types)
{
    inited = true;
    settings_block.copy(&packing_params_block);
    ctx.types = types;
    Block *b;
    b = settings_block.get_block("trunks_params");
    if (b)
        trunks_params.copy(b);
    
    b = settings_block.get_block("branches_params");
    if (b)
        branches_params.copy(b);
    
    b = settings_block.get_block("trees_params");
    if (b)
        trees_params.copy(b);
}

void GrovePacker::prepare_grove_atlas(GrovePacked &grove, int tex_w, int tex_h, bool save_atlases, bool save_png, 
                                      bool alpha_tex_needed)
{
    auto &atl = grove.groveTexturesAtlas;
    atl.leaves_tex_map.clear();
    atl.wood_tex_map.clear();
    atl.atlases_valid = false;
    atl.maps_valid = false;
    if (atl.leavesAtlas)
      delete atl.leavesAtlas;
    if (atl.woodAtlas)
      delete atl.woodAtlas;

    std::vector<Texture> unique_wood_texs;
    std::vector<Texture> unique_leaves_texs;
    std::vector<bool> ids_found = std::vector<bool>(1024, false);
    for (auto &ib : grove.instancedBranches)
    {
        for (int &tid : ib.IDA.type_ids)
        {
            if (!ids_found[tid])
            {
                ids_found[tid] = true;
                auto &type = metainfoManager->get_tree_type(tid);
                int w_tex_n = unique_wood_texs.size();
                for (int i=0;i<unique_wood_texs.size();i++)
                {
                    if (unique_wood_texs[i].texture == type.wood.texture)
                    {
                        w_tex_n = i;
                        break;
                    }
                }
                if (w_tex_n == unique_wood_texs.size())
                    unique_wood_texs.push_back(type.wood);

                int l_tex_n = unique_leaves_texs.size();
                for (int i=0;i<unique_leaves_texs.size();i++)
                {
                    if (unique_leaves_texs[i].texture == type.leaf.texture)
                    {
                        l_tex_n = i;
                        break;
                    }
                }
                if (l_tex_n == unique_leaves_texs.size())
                    unique_leaves_texs.push_back(type.leaf);

                atl.wood_tex_map.emplace(tid, w_tex_n);
                atl.leaves_tex_map.emplace(tid, l_tex_n);
            }
        }
    }
    atl.maps_valid = true;
    atl.atlases_valid = save_atlases;
    atl.woodAtlas = new TextureAtlas(tex_w, tex_h*unique_wood_texs.size(), 1);
    atl.woodAtlas->set_grid(tex_w, tex_h, false);
    atl.leavesAtlas = new TextureAtlas(tex_w, tex_h*unique_leaves_texs.size(), 1);
    atl.leavesAtlas->set_grid(tex_w, tex_h, false);

    PostFx copy = PostFx("copy.fs");

    for (int i=0;i<unique_wood_texs.size();i++)
    {
        int tex_id = atl.woodAtlas->add_tex();
        atl.woodAtlas->target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex",unique_wood_texs[i]);
        copy.render();
        for (auto &p : atl.wood_tex_map)
        {
            if (p.second == i)
                p.second = tex_id;
        }
    }

    if (save_png)
        engine::textureManager->save_png(atl.woodAtlas->tex(0),"wood_atlas");
    if (save_atlases)
        atl.woodAtlas->gen_mipmaps();
    else
    {
        delete atl.woodAtlas;
        atl.woodAtlas = nullptr;
    }

    for (int i=0;i<unique_leaves_texs.size();i++)
    {
        int tex_id = atl.leavesAtlas->add_tex();
        atl.leavesAtlas->target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex",unique_leaves_texs[i]);
        copy.render();
        for (auto &p : atl.leaves_tex_map)
        {
            if (p.second == i)
                p.second = tex_id;
        }
    }

    if (save_png)
        engine::textureManager->save_png(atl.leavesAtlas->tex(0),"leaves_atlas");
    if (alpha_tex_needed)
    {
        PostFx copy_alpha = PostFx("alpha_split_alpha.fs");
        TextureAtlas atl_alpha = TextureAtlas(tex_w, tex_h*unique_leaves_texs.size(), 1);
        atl_alpha.set_grid(tex_w, tex_h, false);

        for (int i=0;i<unique_leaves_texs.size();i++)
        {
            int tex_id = atl_alpha.add_tex();
            atl_alpha.target_slice(tex_id, 0);
            copy_alpha.use();
            copy_alpha.get_shader().texture("tex",unique_leaves_texs[i]);
            copy_alpha.render();
        }
        engine::textureManager->save_png(atl_alpha.tex(0),"leaves_atlas_alpha");
    }
    if (save_atlases)
        atl.leavesAtlas->gen_mipmaps();
    else
    {
        delete atl.leavesAtlas;
        atl.leavesAtlas = nullptr;
    }
}

void GrovePacker::remove_trees(GrovePacked &grove, std::vector<int> &ids)
{
  //tmp structures that connects internal data from grove packer with grove structure
  std::vector<ClusterPackingLayer *> cpls;
  std::vector<int> layer_by_cpl_n;
  for (auto &cpl : packingLayersBranches)
  {
    layer_by_cpl_n.push_back(0);
    cpls.push_back(&cpl);
  }
  for (auto &cpl : packingLayersTrunks)
  {
    layer_by_cpl_n.push_back(1);
    cpls.push_back(&cpl);
  }
  for (auto &cpl : packingLayersTrees)
  {
    layer_by_cpl_n.push_back(2);
    cpls.push_back(&cpl);
  }
  //key is mesh id of first mesh of instanced branch
  //ivec2 is <index of ClusterPackingLayer in cpls>,<index of cluster inside it>
  //we don't modify these arrays until the very end, so it's safe to use that indexing 
  std::map<int, glm::ivec2> cluster_by_packed_branch; 
  int total_cl_cnt = 0;
  int br_cl_cnt = 0;
  for (int i=0;i<cpls.size();i++)
  {
    if (cpls[i]->clusters.size() != cpls[i]->additional_data.size())
          logerr("Clusted packing layer is corrupted. %d clusters, %d additional data", cpls[i]->clusters.size(), cpls[i]->additional_data.size());
    for (int j=0;j<cpls[i]->clusters.size();j++)
    {
      total_cl_cnt++;
      if (cpls[i]->additional_data[j].has_instanced_branch)
      {
        br_cl_cnt++;
        if (!cluster_by_packed_branch.emplace(cpls[i]->additional_data[j].instanced_branch->branches[0], glm::ivec2(i,j)).second)
        {
          logerr("Clusted %d %d connected to the already existed branch",i,j);
        }
      }
      else
      {
        //it is not a problem, if cluster represent impostor, it does not have an instanced branch
        //logerr("Clusted %d %d does not have branch",i,j);
      }
    }
  }
  if (br_cl_cnt != grove.instancedBranches.size())
    logerr("we have %d clusters, %d with branches and %d instanced branches",total_cl_cnt, br_cl_cnt, grove.instancedBranches.size());
  for (auto it=grove.instancedBranches.begin(); it != grove.instancedBranches.end(); it++)
  {
    auto mit = cluster_by_packed_branch.find(it->branches[0]);
    if (mit==cluster_by_packed_branch.end())
      logerr("instanced branch is not presented in any cluster");
  }
  //bitsets representing what to delete
  std::vector<bool> models_to_delete(grove.instancedBranches.size(), false);
  std::vector<std::vector<bool>> instances_to_delete(grove.instancedBranches.size());
  std::vector<int> instances_detele_cnt(grove.instancedBranches.size(), 0);
  {
  int i=0;
  for (auto &m : grove.instancedBranches)
  {
    instances_to_delete[i] = std::vector<bool>(m.IDA.centers_par.size(), false);
    i++;
  }
  }

  //filling these bitsets
  for (int id : ids)
  {
    auto it = grove.trees_by_global_id.find(id);
    if (it == grove.trees_by_global_id.end())
      continue;
    
    std::vector<CompressedTree::Node> nodes = grove.compressedTrees[it->second].LOD_roots;
    while (!nodes.empty())
    {
      std::vector<CompressedTree::Node> new_nodes;
      for (auto &node : nodes)
      {
        if (node.type == CompressedTree::MODEL)
        {
          if (instances_to_delete[node.model_num][node.instance_num] == false)
            instances_detele_cnt[node.model_num]++;
          instances_to_delete[node.model_num][node.instance_num] = true;
        }
        new_nodes.insert(new_nodes.end(), node.children.begin(), node.children.end());
      }
      nodes = new_nodes;
    }
  }

  {
  int i=0;
  for (auto &m : grove.instancedBranches)
  {
    if (instances_detele_cnt[i] == m.IDA.centers_par.size())
      models_to_delete[i] = true;
    i++;
  }
  }

  int removed_ib = 0;
  //remove all needed structures in efficient way
  {
  int i=0;
  auto it = grove.instancedBranches.begin();
  while (it != grove.instancedBranches.end())
  {
    glm::ivec2 l_cl = cluster_by_packed_branch[it->branches[0]];
    auto &ACDA = cpls[l_cl.x]->clusters[l_cl.y].ACDA;
    // we need to remove clustering data using right clusterizer
    //(because there is different data inside depends of clusterizer)
    // don't worry - creating of a clusterizer is fast, as it does not contain heavy data
    Clusterizer2 *cl = nullptr;
    if (layer_by_cpl_n[l_cl.x] == 0)
      cl = get_clusterizer(branches_params, cStrategy);
    else if (layer_by_cpl_n[l_cl.x] == 1)
      cl = get_clusterizer(trunks_params, cStrategy);
    else
      cl = get_clusterizer(trees_params, cStrategy);

    //if we removed all instances of this branch, we can delete it
    if (models_to_delete[i])
    {
      //we do not remove clusters immediately, because it will spoil cluster_by_packed_branch structure
      cpls[l_cl.x]->clusters[l_cl.y].base = nullptr;//mark as invalid
      removed_ib++;
      //then remove all meshes of this branch
      for (int id : it->branches)
        grove.instancedCatalogue.remove(id);
      
      //then remove the branch itself
      it = grove.instancedBranches.erase(it);

      //we still need to clear all clustering data from ACDA
      for (BranchClusteringData *cd : ACDA.clustering_data)
      {
        if (cd)
          cl->clear_branch_data(cd, &ctx);
      }
    }
    else
    {
      //this branch remains, but we probably need to delete some of it's instances
      for (int j=it->IDA.centers_par.size()-1;j>=0;j--)
      {
        if (instances_to_delete[i][j])
        {
          //remove instances from IDA
          it->IDA.centers_par.erase(it->IDA.centers_par.begin()+j);
          it->IDA.centers_self.erase(it->IDA.centers_self.begin()+j);
          it->IDA.transforms.erase(it->IDA.transforms.begin()+j);
          it->IDA.tree_ids.erase(it->IDA.tree_ids.begin()+j);
          it->IDA.type_ids.erase(it->IDA.type_ids.begin()+j);

          //remove instances from ACDA
          if (ACDA.clustering_data[j])
          { 
            cl->clear_branch_data(ACDA.clustering_data[j], &ctx);
          }
          ACDA.clustering_data.erase(ACDA.clustering_data.begin()+j);
        }
      }
      it++;
    }
    delete cl;
    i++;
  }
  }

  int removed_cl = 0;
  //remove unnecessary clusters
  for (auto *cpl : cpls)
  {
    for (int i=cpl->additional_data.size()-1;i>=0;i--)
    {
      if (cpl->clusters[i].base == nullptr)
      {
        cpl->clusters.erase(cpl->clusters.begin()+i);
        cpl->additional_data.erase(cpl->additional_data.begin()+i);
        removed_cl++;
      }
    }
  }
  if (removed_ib != removed_cl)
    logerr("removed %d instanced branches and %d clusters", removed_ib, removed_cl);
  //recalculate tree structures as they contain positions of instances that may change
  recreate_compressed_trees(grove);
}

void GrovePacker::clear()
{
  originalBranches.branches.clear();
  originalLeaves.leaves.clear();

  std::vector<ClusterPackingLayer *> cpls;
  std::vector<int> layer_by_cpl_n;
  for (auto &cpl : packingLayersBranches)
  {
    layer_by_cpl_n.push_back(0);
    cpls.push_back(&cpl);
  }
  for (auto &cpl : packingLayersTrunks)
  {
    layer_by_cpl_n.push_back(1);
    cpls.push_back(&cpl);
  }
  for (auto &cpl : packingLayersTrees)
  {
    layer_by_cpl_n.push_back(2);
    cpls.push_back(&cpl);
  }

  for (int i=0;i<layer_by_cpl_n.size();i++)
  {
    int l = layer_by_cpl_n[i];
    Clusterizer2 *cl = nullptr;
    if (layer_by_cpl_n[l] == 0)
      cl = get_clusterizer(branches_params, cStrategy);
    else if (layer_by_cpl_n[l] == 1)
      cl = get_clusterizer(trunks_params, cStrategy);
    else
      cl = get_clusterizer(trees_params, cStrategy);
    for (auto &cluster : cpls[i]->clusters)
    {
      for (BranchClusteringData *cd : cluster.ACDA.clustering_data)
      {
        if (cd)
          cl->clear_branch_data(cd, &ctx);
      }
    }
    delete cl;
  }  

  packingLayersBranches.clear();
  packingLayersTrees.clear();
  packingLayersTrunks.clear();
  ctx.clear();

  ib_id_counter = 1;
}