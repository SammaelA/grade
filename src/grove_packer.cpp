#include "tinyEngine/utility.h"
#include "clustering/clustering.h"
#include "texture_manager.h"
#include "visualizer.h"
#include "distribution.h"
#include <math.h>
#include <algorithm>
#include "body.h"
#include <chrono>
#include "tinyEngine/save_utils/saver.h"
#include "impostor.h"
#include "terrain.h"
#include "field_2d.h"
#include "grove_generation_utils.h"
#include "synthetic_trees_generator.h"
#include "grove_packer.h"
#include "clustering/default_clustering_params.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

GrovePacker::GrovePacker(bool shared_ctx)
{
    shared_context = shared_ctx;
}

void pack_branch_recursively(::Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, BranchStructure &b_struct, int lvl_from, int lvl_to);
std::list<InstancedBranch>::iterator pack_cluster(ClusterData &cluster, GrovePacked &grove,
                                                  std::vector<BranchStructure> &instanced_structures, int lvl_from, int lvl_to)
{
    instanced_structures.push_back(BranchStructure());
    grove.instancedBranches.push_back(InstancedBranch());
    grove.instancedBranches.back().IDA = cluster.IDA;
    std::vector<unsigned> &ids = grove.instancedBranches.back().branches;
    ::Branch *base = cluster.base;
    pack_branch_recursively(base, grove, ids, instanced_structures.back(), lvl_from, lvl_to);
    if (instanced_structures.back().childBranches.size() == 1)
    {
        BranchStructure bs = instanced_structures.back().childBranches[0];
        instanced_structures.back() = bs;
    }
    for (int i = 0; i < cluster.ACDA.originals.size(); i++) //leave marks on branch to construct tree structure in future
    {
        if (cluster.ACDA.originals[i])
        {
            cluster.ACDA.originals[i]->mark_A = instanced_structures.size() - 1; //cluster id
            cluster.ACDA.originals[i]->mark_B = -i - 100;                        //i is number in cluster
        }
    }
    grove.instancedBranches.back().bbox = BillboardCloudRaw::get_minimal_bbox(base);

    std::list<InstancedBranch>::iterator it = grove.instancedBranches.end();
    it--;
    return it;
}

void pack_structure(::Branch *rt, GrovePacked &grove, BranchStructure &str, std::vector<BranchStructure> &instanced_structures)
{
    return;
    /*
    str.pos = rt->mark_B;
    for (::Joint &j : rt->joints)
    {
        for (::Branch *br : j.childBranches)
        {
            if (br->mark_B < 0)//it is a mark made by pack_cluster
            {
                unsigned transform_n = - (br->mark_B + 100);
                unsigned instance_n = br->mark_A;
                BranchStructure bs = instanced_structures[instance_n];
                glm::mat4 tr = grove.instancedBranches[instance_n].IDA.transforms[transform_n];
                str.childBranchesInstanced.push_back(std::pair<glm::mat4,BranchStructure>(tr,bs));
            }
            else
            {
                BranchStructure bs;
                pack_structure(br,grove,bs,instanced_structures);
                str.childBranches.push_back(bs);
            }
            
        }
    }*/
}
void pack_branch_recursively(::Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, BranchStructure &b_struct, int lvl_from, int lvl_to)
{
    if (b->level > lvl_to)
        return;
    if (b->level >= lvl_from)
    {
        PackedBranch pb;
        b->pack(pb);
        ids.push_back(grove.instancedCatalogue.add(pb, b->level));
        //b_struct.childBranches.push_back(BranchStructure(ids.back()));
    }
    for (::Joint &j : b->joints)
    {
        for (::Branch *br : j.childBranches)
            pack_branch_recursively(br, grove, ids, b_struct, lvl_from, lvl_to);
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
bool is_valid_tree(::Tree &t)
{
    if (!t.valid)
    {
        logerr("tree %u was marked as not valid by a generator. Grove packer will ignore it", t.id);
        return false;
    }
    if (!(t.valid && t.leaves && t.branchHeaps.size() >= 2))
    {
        logerr("tree %u do not have some of essential data structures. %d Grove packer will ignore it", t.id,
               t.branchHeaps.size());
        return false;
    }
    if (t.branchHeaps[0]->branches.empty())
    {
        logerr("tree %u do not have trunk. Grove packer will ignore it", t.id);
        return false;
    }
    if (t.branchHeaps[1]->branches.empty())
    {
        logerr("tree %u do not have branches instead of trunk. Grove packer will ignore it", t.id);
        return false;
    }
    return true;
}

void GrovePacker::pack_layer(Block &settings, GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                             std::vector<ClusterPackingLayer> &packingLayers, LightVoxelsCube *post_voxels,
                             int layer_from, int layer_to, bool models, bool bill, bool imp,
                             bool visualize_clusters)
{


    std::string c_strategy_name = "merge";
    c_strategy_name = settings.get_string("clustering_strategy",c_strategy_name);

    if (c_strategy_name == "recreate")
    {
        cStrategy = ClusteringStrategy::Recreate;
    }
    else 
    {
        cStrategy = ClusteringStrategy::Merge;
    }

    bool need_something = ggd.task & (GenerationTask::SYNTS | GenerationTask::CLUSTERIZE) ||
                          ((ggd.task & GenerationTask::MODELS) && models) ||
                          ((ggd.task & GenerationTask::BILLBOARDS) && bill) ||
                          ((ggd.task & GenerationTask::IMPOSTORS) && imp);
    if (!need_something)
    {
        return;
    }

    int count = ggd.trees_count;
    int clusters_before = packingLayers[0].clusters.size();
    Clusterizer2 *cl = new Clusterizer2(cStrategy);

    ctx->light = post_voxels;
    ctx->types = &(ggd.types);
    cl->prepare(settings);

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
    if ((ggd.task & (GenerationTask::CLUSTERIZE)))
    {
        std::vector<ClusterData> clusters_base;
        //create base clusters one element in each
        cl->get_base_clusters(settings, trees_external,count, layer_from, clusters_base, ctx);
        //clusterize this base clusters and put result at the end of cluster list
        //on zero level
        cl->clusterize(settings, clusters_base, packingLayers[clustering_base_level].clusters, ctx, 
                       save_clusterizer, visualize_clusters);
        if (save_clusterizer)
        {
            saved_clustering_data.push_back(cl->get_full_data());
        }
    }
    else
    {
        //just add base clusters to zero level cluster list
        cl->get_base_clusters(settings, trees_external,count, layer_from, packingLayers[clustering_base_level].clusters, ctx);
    }

    //we do not modify previously existed clusters on this step. Only add some new
    for (int i = clusters_before; i < packingLayers[clustering_base_level].clusters.size(); i++)
    {
        packingLayers[clustering_base_level].additional_data.emplace_back();
        new_clusters.push_back(ClusterInfo{clustering_base_level, i});
    }

    if (ggd.task & (GenerationTask::CLUSTERIZE))
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
                cl->clusterize(settings, packingLayers[i].clusters, packingLayers[i + 1].clusters, ctx,
                               save_clusterizer, visualize_clusters);
                int new_size = packingLayers[i + 1].clusters.size();

                for (int j = 0; j < prev_size; j++)
                {
                    old_cl_poses.emplace(packingLayers[i].clusters[j].id, j);
                    debugl(6, "old cluster %d %d\n", j, (int)(packingLayers[i].clusters[j].id));
                }
                for (int j = 0; j < new_size; j++)
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
    
        if (models && (ggd.task & (GenerationTask::MODELS)))
        {
            std::vector<BranchStructure> instanced_structures;

            auto it = pack_cluster(packingLayers[info.layer].clusters[info.pos], grove, instanced_structures, layer_from, layer_to);
            packingLayers[info.layer].additional_data[info.pos].instanced_branch = it;
            packingLayers[info.layer].additional_data[info.pos].is_presented = true;
        }
        if (bill && (ggd.task & (GenerationTask::BILLBOARDS)))
        {
            BillboardCloudRaw cloud;
            cloud.prepare(ULTRALOW, layer_from, packingLayers[info.layer].clusters[info.pos], ggd.types,
                          &grove.clouds[2], packingLayers[info.layer].additional_data[info.pos].large_billboards);
            cloud.prepare(ULTRALOW, layer_from + 1, packingLayers[info.layer].clusters[info.pos], ggd.types,
                          &grove.clouds[3], packingLayers[info.layer].additional_data[info.pos].small_billboards);
        }
        if (imp && (ggd.task & (GenerationTask::IMPOSTOR_FULL_GROVE)))
        {
            debugl(6, "full grove impostors are temporary unavailable\n");
        }
        if (imp && (ggd.task & (GenerationTask::IMPOSTORS)))
        {
            ImpostorBaker ib;
            ib.prepare(ggd.impostor_quality, layer_from, packingLayers[info.layer].clusters[info.pos], ggd.types,
                       &(grove.impostors[1]), packingLayers[info.layer].additional_data[info.pos].impostors);
        }

        debugl(6, "added cluster %d\n", packingLayers[info.layer].clusters[info.pos].id);
    }
    for (auto &info : expanded_clusters)
    {       
        if (layer_from == 0 && layer_to == 0)
            recalculate_nodes(packingLayers[info.to.layer].clusters[info.to.pos]); 
        else if (layer_from == 1)
            transform_by_nodes(packingLayers[info.to.layer].clusters[info.to.pos]);
        
        if (models && (ggd.task & (GenerationTask::MODELS)))
        {
            auto &it = packingLayers[info.from.layer].additional_data[info.from.pos].instanced_branch;
            it->IDA = packingLayers[info.to.layer].clusters[info.to.pos].IDA;
            packingLayers[info.to.layer].additional_data[info.to.pos].instanced_branch = it;
            packingLayers[info.to.layer].additional_data[info.to.pos].is_presented = true;
        }
        if (bill && (ggd.task & (GenerationTask::BILLBOARDS)))
        {
            BillboardCloudRaw cloud;
            auto &from_a = packingLayers[info.from.layer].additional_data[info.from.pos];
            auto &to_a = packingLayers[info.to.layer].additional_data[info.to.pos];
            to_a.large_billboards = from_a.large_billboards;
            to_a.small_billboards = from_a.small_billboards;

            cloud.extend(ULTRALOW, layer_from, packingLayers[info.to.layer].clusters[info.to.pos], ggd.types,
                         &grove.clouds[2], to_a.large_billboards);
            cloud.extend(ULTRALOW, layer_from + 1, packingLayers[info.to.layer].clusters[info.to.pos], ggd.types,
                         &grove.clouds[3], to_a.small_billboards);

            packingLayers[info.to.layer].additional_data[info.to.pos].is_presented = true;
        }
        if (imp && (ggd.task & (GenerationTask::IMPOSTORS)))
        {
            auto &from_a = packingLayers[info.from.layer].additional_data[info.from.pos];
            auto &to_a = packingLayers[info.to.layer].additional_data[info.to.pos];

            to_a.impostors = from_a.impostors;
            to_a.impostors->IDA = packingLayers[info.to.layer].clusters[info.to.pos].IDA;
        }
        debugl(6, "replaced cluster %d with %d\n", packingLayers[info.from.layer].clusters[info.from.pos].id,
               packingLayers[info.to.layer].clusters[info.to.pos].id);
    }
    for (auto &info : removed_clusters)
    {
        bool need_clear_data = true;
        if (need_clear_data)
        {
            if (models && (ggd.task & (GenerationTask::MODELS)))
            {
                auto &it = packingLayers[info.layer].additional_data[info.pos].instanced_branch;
                for (unsigned &id : it->branches)
                    grove.instancedCatalogue.remove(id);
                grove.instancedBranches.erase(it);
            }
            if (bill && (ggd.task & (GenerationTask::BILLBOARDS)))
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
            if (imp && (ggd.task & (GenerationTask::IMPOSTORS)))
            {
                auto &it = packingLayers[info.layer].additional_data[info.pos].impostors;
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
                glm::vec3 dist = node.position - pos;
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

void GrovePacker::add_trees_to_grove(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                                     bool visualize_clusters, bool save_cluster_data)
{
    int max_trees_in_patch = settings_block.get_int("max_trees_in_patch",1000);
    if (ggd.trees_count <= max_trees_in_patch)
    {
        add_trees_to_grove_internal(ggd, grove, trees_external, h, visualize_clusters, save_cluster_data);
    }
    else
    {
        int start_pos = 0;
        int trees_count = ggd.trees_count;
        ggd.trees_count = max_trees_in_patch;
        while (start_pos < trees_count)
        {
            add_trees_to_grove_internal(ggd, grove, trees_external + start_pos, h, visualize_clusters, save_cluster_data);
            start_pos+=max_trees_in_patch;
        }
        ggd.trees_count = trees_count;
    }
}

void GrovePacker::add_trees_to_grove_prepare_dataset(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, 
                                                     Heightmap *h, std::string &save_path)

{
    add_trees_to_grove_internal(ggd, grove, trees_external, h, false, true);

    int cnt = 0;
    TextureAtlasRawData raw_atlas = TextureAtlasRawData(ctx->self_impostors_data->atlas);
    unsigned char *sl_data = safe_new<unsigned char>(2*raw_atlas.get_slice_size(0), "sl_data");
    memset(sl_data,0,2*raw_atlas.get_slice_size(0));
    unsigned char *tmp_data = safe_new<unsigned char>(raw_atlas.get_slice_size(0), "tmp_data");
    int ww = 0, hh = 0;
    bool info_files = true;
    
    std::string database_file;
    std::string test_file;
    std::string train_file;
    
    float train_part = 0.9;
    float test_part = 0.05;
    
    std::string database_name = "database.txt";
    std::string test_name = "test.txt";
    std::string train_name = "train.txt";
    
    std::string folder_name = "images";

    int sl = 0;
    int train_elems = 0;
    int test_elems = 0;

    std::string dir_path;

    //find largest branch to rescale dataset images properly
    glm::vec3 max_sizes = glm::vec3(0,0,0);
    for (int i = 0; i< packingLayersBranches.size();i++)
    {
        for (auto &c : packingLayersBranches[i].clusters)
        {
            for (auto *cd : c.ACDA.clustering_data)
            {
                auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                if (imp_cd)
                {
                    max_sizes = max(max_sizes, imp_cd->min_bbox.sizes);
                }
            }
        }
    }
    float q = 0.85;
    float max_size = q*MAX(max_sizes.x,MAX(max_sizes.y,max_sizes.z));
    float need_rescale = true;
    try
    {
        int clusters_count = 0;
        if (info_files)
        {
            dir_path = save_path + "/"+folder_name;
            boost::filesystem::create_directory(dir_path);
            boost::filesystem::permissions(dir_path, boost::filesystem::perms::all_all); 

            for (int i = 0; i< packingLayersBranches.size();i++)
            {
                clusters_count += packingLayersBranches[i].clusters.size();
            }         
        }
        for (int i = 0; i< packingLayersBranches.size();i++)
        {
            for (auto &c : packingLayersBranches[i].clusters)
            {
                std::string cluster_labels;
                if (info_files)
                {
                    for (int j=0;j<clusters_count;j++)
                    {
                        cluster_labels += ((j == cnt) ? " 1" : " 0");
                    }
                }
                else
                {
                    dir_path = save_path + "/" + std::to_string(cnt);
                    boost::filesystem::create_directory(dir_path);
                    boost::filesystem::permissions(dir_path, boost::filesystem::perms::all_all);
                    sl = 0;
                }
                for (auto *cd : c.ACDA.clustering_data)
                {
                    auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                    if (imp_cd)
                    {
                        for (auto &bill : imp_cd->self_impostor->slices)
                        {
                            std::string file_path = dir_path + "/" + std::to_string(sl)+".bmp";
                            raw_atlas.get_slice(bill.id, sl_data, &ww, &hh);

                            if (need_rescale)
                            {
                                float sz = MAX(imp_cd->min_bbox.sizes.x,MAX(imp_cd->min_bbox.sizes.y,imp_cd->min_bbox.sizes.z));
                                float scale = max_size/sz;
                                if (true)
                                {
                                    for (int y=0;y<hh;y++)
                                    {
                                        float y_src = scale*y;
                                        int y0 = y_src;
                                        float qy = y_src - y0;
                                        for (int x=0;x<ww;x++)
                                        {
                                            float x_src = scale*x;
                                            if (x_src > ww + 1 || y_src > hh + 1)
                                            {
                                                for (int ch = 0;ch < 3;ch++)
                                                    tmp_data[4*(y*ww + x) + ch] = 0;
                                            }
                                            else
                                            {
                                                int x0 = x_src;
                                                float qx = x_src - x0;
                                                for (int ch = 0;ch < 3;ch++)
                                                {
                                                    tmp_data[4*(y*ww + x) + ch] = 
                                                             (1-qy)*((1-qx)*sl_data[4*(y0*ww + x0) + ch] + 
                                                                      qx*sl_data[4*(y0*ww + x0 + 1) + ch]) +
                                                                qy*((1-qx)*sl_data[4*((y0+1)*ww + x0) + ch] + 
                                                                      qx*sl_data[4*((y0+1)*ww + x0 + 1) + ch]);
                                                    tmp_data[4*(y*ww + x) + ch] = sl_data[4*(y0*ww + x0) + ch];
                                                }
                                            }
                                            tmp_data[4*(y*ww + x) + 3] = 255;
                                        }
                                    }
                                    textureManager.save_bmp_raw_directly(tmp_data, ww, hh, 4, file_path);
                                }
                                else
                                {
                                    textureManager.save_bmp_raw_directly(sl_data, ww, hh, 4, file_path);
                                }
                            }
                            else
                                textureManager.save_bmp_raw_directly(sl_data, ww, hh, 4, file_path);
                            if (info_files)
                            {
                                //add a record about images to database and (maybe) test or train lists;
                                std::string record = folder_name + "/" + std::to_string(sl)+".bmp" + cluster_labels + "\n";
                                database_file += record;
                                if (urand() < train_part)
                                {
                                    train_file += record;
                                    train_elems++;
                                }
                                if (urand() < test_part)
                                {
                                    test_file += record;
                                    test_elems++;
                                }
                            }
                            sl++;
                        }
                    }
                    else
                    {
                        logerr("error - trying to save to dataset branch without impostor. Check clustering settings");
                    }
                }
                cnt++;
                if (info_files)
                {
                    std::ofstream database_ofs;
                    database_ofs.open(save_path+"/database.txt");
                    database_ofs << database_file;
                    database_ofs.close();

                    std::ofstream test_ofs;
                    test_ofs.open(save_path+"/test.txt");
                    test_ofs << test_file;
                    test_ofs.close();

                    std::ofstream train_ofs;
                    train_ofs.open(save_path+"/train.txt");
                    train_ofs << train_file;
                    train_ofs.close();

                    std::ofstream info_ofs;
                    info_ofs.open(save_path+"/info.txt");
                    info_ofs << std::string(std::to_string(clusters_count)+" "+std::to_string(sl)+" "+
                                            std::to_string(train_elems)+" "+std::to_string(test_elems));
                    info_ofs.close();

                }
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    safe_delete<unsigned char>(sl_data, "sl_data");
    safe_delete<unsigned char>(tmp_data, "tmp_data");
    raw_atlas.clear();
}

void GrovePacker::add_trees_to_grove_internal(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                                              bool visualize_clusters, bool save_cluster_data)
{
    if (!inited)
        init();
    
    for (int i = grove.impostors.size(); i < 2; i++)
    {
        grove.impostors.push_back(ImpostorsData());
        grove.impostors.back().valid = true;
    }
    for (int i = grove.clouds.size(); i < 4; i++)
        grove.clouds.push_back(BillboardCloudData());

    save_clusterizer = save_cluster_data;
    grove.center = glm::vec3(0, 0, 0);
    grove.ggd_name = ggd.name;
    int synts = ggd.synts_count;
    int count = ggd.trees_count;
    LightVoxelsCube *post_voxels = nullptr;
    Seeder *post_seeder = nullptr;
    GroveGenerationData &curGgd = ggd;

    int valid_trees_cnt = 0;
    for (int i = 0; i < count; i++)
    {
        bool valid = is_valid_tree(trees_external[i]);
        trees_external[i].valid = valid;
        valid_trees_cnt += valid;
    }
    if (count == 0)
    {
        logerr("Grove %s is empty. ", ggd.name.c_str());
        return;
    }
    if (valid_trees_cnt == 0)
    {
        logerr("Grove %s has %d tree(s) but none of them are valid.", ggd.name.c_str(), count);
        return;
    }

    debugl(1, "Packing grove %s with %d/%d valid trees\n", ggd.name.c_str(), valid_trees_cnt, count);

    {
        float r = sqrt(count);
        glm::vec3 vox_center = glm::vec3(0, 100, 0) + curGgd.pos;
        glm::vec3 vox_size = curGgd.size;
        auto &type = ggd.types[0];
        TreeStructureParameters base_params;
        ParameterSetWrapper params = ParameterSetWrapper(base_params, base_params.max_depth() + 1);
        params.set_state(params().max_depth() - 1);
        float single_voxel_size = params().seg_len_mult() / params().light_precision();

        post_voxels = new LightVoxelsCube(vox_center, vox_size, params().seg_len_mult(), params().light_precision());
        post_seeder = new Seeder(ggd, 10, h);

        post_voxels->add_heightmap(*h);
        for (int i = 0; i < curGgd.obstacles.size(); i++)
        {
            post_voxels->add_body(curGgd.obstacles[i]);
            post_seeder->add_body(curGgd.obstacles[i]);
        }
        post_seeder->recalcuate_shadows(trees_external, count);
        add_occluder(post_voxels, trees_external, count);
    }

    for (int i = 0; i < count; i++)
    {
        auto it = trees_nodes.emplace(trees_external[i].id,std::vector<Node>{});
        for (Joint &j : trees_external[i].root->joints)
        {
            it.first->second.push_back({j.pos});
        }
    }


    current_clustering_step = ClusteringStep::TRUNKS;
    pack_layer(*trunks_params, ggd, grove, trees_external, h, packingLayersTrunks, post_voxels,
               0, 0, true, false, false, false);

    current_clustering_step = ClusteringStep::BRANCHES;
    pack_layer(*branches_params, ggd, grove, trees_external, h, packingLayersBranches, post_voxels,
               1, 1000, true, true, false, visualize_clusters);

    current_clustering_step = ClusteringStep::TREES;
    pack_layer(*trees_params, ggd, grove, trees_external, h, packingLayersTrees, post_voxels,
               0, 1000, false, false, true, false);

    //transform_all_according_to_root(grove);
    
    originalBranches.clear_removed();
    originalLeaves.clear_removed();

    delete (post_voxels);
    delete (post_seeder);
}
void GrovePacker::base_init()
{
    if (shared_context)
        ctx = new ClusteringContext();
    else    
        ctx = &self_ctx;
    Block *b;
    b = settings_block.get_block("trunks_params");
    if (b)
        trunks_params = b;
    
    b = settings_block.get_block("branches_params");
    if (b)
        branches_params = b;
    
    b = settings_block.get_block("trees_params");
    if (b)
        trees_params = b;
}
void GrovePacker::init(Block &packing_params_block)
{
    inited = true;
    settings_block = packing_params_block;
    base_init();
}
void GrovePacker::init()
{
    inited = true;
    BlkManager man;
    man.load_block_from_file("settings.blk",settings_block);  
    base_init(); 
}
