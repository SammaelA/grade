#include "tinyEngine/utility.h"
#include "branch_clusterization.h"
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

void transform_according_to_root(ClusterData &cluster)
{
    InstanceDataArrays IDA = cluster.IDA;
    ::Branch *base = cluster.base;
    if (cluster.ACDA.originals.size() != IDA.transforms.size())
    {
        logerr("Transform according to root failed: corrupted clusters data.");
        return;
    }
    for (int i = 0; i < cluster.ACDA.originals.size(); i++)
    {
        if (base->joints.empty() || cluster.ACDA.originals[i]->joints.empty())
            continue;
        std::vector<glm::vec3> replaced_joints;
        for (::Joint &j : base->joints)
            replaced_joints.push_back(glm::vec3(IDA.transforms[i] * glm::vec4(j.pos, 1)));
        for (::Joint &j : cluster.ACDA.originals[i]->joints)
        {
            float d_min = 1e9;
            glm::vec3 n_min = glm::vec3(0);
            for (glm::vec3 &np : replaced_joints)
            {
                float d = glm::length(np - j.pos);
                if (d < d_min)
                {
                    d_min = d;
                    n_min = np;
                }
            }
            glm::vec3 sh = n_min - j.pos;
            for (::Branch *ch_b : j.childBranches)
            {
                glm::mat4 tr = glm::translate(glm::mat4(1.0f), sh);
                ch_b->transform(tr);
            }
        }
    }
}

void GrovePacker::transform_all_according_to_root(GrovePacked &grove)
{
    const int needed_level = 1;
    struct TreeInstance
    {
        std::vector<std::vector<PackedJoint> *> joints;
        glm::mat4 transform;
    };
    std::map<int, TreeInstance> all_trees;

    for (auto &layer : packingLayersTrunks)
    {
        for (int i = 0; i < layer.clusters.size(); i++)
        {
            auto &cl = layer.clusters[i];
            auto &ids = layer.additional_data[i].instanced_branch->branches;
            std::vector<std::vector<PackedJoint> *> joints;

            for (uint br_id : ids)
            {
                auto &f = grove.instancedCatalogue.get(br_id);
                joints.push_back(&(f.joints));
            }

            for (int j = 0; j < cl.IDA.tree_ids.size(); j++)
            {
                TreeInstance in;
                in.joints = joints;
                in.transform = cl.IDA.transforms[j];
                all_trees.emplace(cl.IDA.tree_ids[j], in);
            }
        }
    }

    for (auto &in_br : grove.instancedBranches)
    {
        PackedBranch &base = grove.instancedCatalogue.get(in_br.branches[0]);
        if (base.level != needed_level)
            continue;
        for (int j = 0; j < in_br.IDA.tree_ids.size(); j++)
        {
            glm::vec3 pos = glm::vec3(in_br.IDA.transforms[j] * glm::vec4(base.joints[0].pos, 1));
            auto tree_it = all_trees.find(in_br.IDA.tree_ids[j]);
            if (tree_it == all_trees.end())
            {
                logerr("found branch from tree with unknow id = %d", in_br.IDA.tree_ids[j]);
                continue;
            }
            glm::vec3 nearest_pos = pos + glm::vec3(0, 1000, 0);
            float nearestDistSq = 1000;
            for (std::vector<PackedJoint> *v1 : tree_it->second.joints)
            {
                for (PackedJoint &j : *v1)
                {
                    glm::vec3 jpos = glm::vec3(tree_it->second.transform * glm::vec4(j.pos, 1));
                    float distSq = glm::dot(jpos - pos, jpos - pos);
                    if (distSq < nearestDistSq)
                    {
                        nearestDistSq = distSq;
                        nearest_pos = jpos;
                    }
                }
            }
            glm::vec3 shift = nearest_pos - pos;
            in_br.IDA.transforms[j] = glm::translate(glm::mat4(1.0f), shift) * in_br.IDA.transforms[j];
            //logerr("translated branch %f %f %f", shift.x, shift.y, shift.z);
        }
    }

    for (auto &bc : grove.clouds)
    {
        for (auto &in_br : bc.billboards)
        {
            for (int j = 0; j < in_br.IDA.tree_ids.size(); j++)
            {
                glm::vec3 pos = glm::vec3(in_br.IDA.transforms[j] * glm::vec4(in_br.base_position, 1));
                auto tree_it = all_trees.find(in_br.IDA.tree_ids[j]);
                if (tree_it == all_trees.end())
                {
                    logerr("found branch from tree with unknow id = %d", in_br.IDA.tree_ids[j]);
                    continue;
                }
                glm::vec3 nearest_pos = pos + glm::vec3(0, 1000, 0);
                float nearestDistSq = 1000;
                for (std::vector<PackedJoint> *v1 : tree_it->second.joints)
                {
                    for (PackedJoint &j : *v1)
                    {
                        glm::vec3 jpos = glm::vec3(tree_it->second.transform * glm::vec4(j.pos, 1));
                        float distSq = glm::dot(jpos - pos, jpos - pos);
                        if (distSq < nearestDistSq)
                        {
                            nearestDistSq = distSq;
                            nearest_pos = jpos;
                        }
                    }
                }
                glm::vec3 shift = nearest_pos - pos;
                in_br.IDA.transforms[j] = glm::translate(glm::mat4(1.0f), shift) * in_br.IDA.transforms[j];
                //logerr("translated bill %f %f %f", shift.x, shift.y, shift.z);
            }
        }
    }
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
        cluster.ACDA.originals[i]->mark_A = instanced_structures.size() - 1; //cluster id
        cluster.ACDA.originals[i]->mark_B = -i - 100;                        //i is number in cluster
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
        //pack_branch_recursively(br, grove, ids, b_struct.childBranches.back(), lvl_from, lvl_to);
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

void GrovePacker::pack_layer(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h,
                std::vector<ClusterPackingLayer> &packingLayers, LightVoxelsCube *post_voxels,
                ClusterizationParams cl_p, int layer_from, int layer_to, bool models, bool bill, bool imp)
{
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
    Clusterizer tr_cl;
    tr_cl.set_light(post_voxels);

    if (ggd.task & (GenerationTask::CLUSTERIZE))
    {
        std::vector<ClusterData> clusters_base;
        tr_cl.get_base_clusters(trees_external, count, layer_from, clusters_base);
        tr_cl.clusterize(cl_p, clusters_base, packingLayers[0].clusters);
    }
    else
    {
        tr_cl.get_base_clusters(trees_external, count, layer_from, packingLayers[0].clusters);
    }

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
    for (int i = clusters_before; i < packingLayers[0].clusters.size(); i++)
    {
        packingLayers[0].additional_data.emplace_back();
        new_clusters.push_back(ClusterInfo{0, i});
    }

    if (ggd.task & (GenerationTask::CLUSTERIZE))
    {
        int max_clusters_in_layer = settings_block.get_int("max_clusters_in_layer",1000);
        int max_layer = settings_block.get_int("recursive_clustering_layers",1);

        for (int i = 0; i < packingLayers.size(); i++)
        {
            int prev_size = packingLayers[i].clusters.size();
            if (i < max_layer && prev_size > max_clusters_in_layer)
            {
                //TODO trunks recursive clusterization - transform according to root
                cleared_layers.push_back(i);
                if (i + 1 == packingLayers.size())
                {
                    packingLayers.push_back(ClusterPackingLayer());
                }
                std::map<long, int> old_cl_poses;
                BitVector remains;
                remains.resize(prev_size, false);
                tr_cl.clusterize(cl_p, packingLayers[i].clusters, packingLayers[i + 1].clusters);
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

        first = false;
        debugl(6, "added cluster %d\n", packingLayers[info.layer].clusters[info.pos].id);
    }
    for (auto &info : expanded_clusters)
    {
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
        debugl(6, "removed cluster %d\n", packingLayers[info.layer].clusters[info.pos].id);
    }
    for (auto i : cleared_layers)
    {
        packingLayers[i].additional_data.clear();
        packingLayers[i].clusters.clear();
        debugl(6, "level cleared %d\n", i);
    }
}

void GrovePacker::add_trees_to_grove(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h)
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

    ClusterizationParams tr_cp;
    tr_cp.weights = std::vector<float>{1, 0, 0, 0.0, 0.0};
    tr_cp.ignore_structure_level = 1;
    tr_cp.delta = 0.1;
    tr_cp.light_importance = 0;
    tr_cp.different_types_tolerance = true;
    tr_cp.r_weights = std::vector<float>{0.4, 0, 0, 0.0, 0.0};
    tr_cp.max_individual_dist = 0.0;
    tr_cp.bwd_rotations = 4;
    //tr_cp.load_from_block(settings_block.get_block("trunk_clusterization_params"));

    ClusterizationParams br_cp;
    br_cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
    br_cp.ignore_structure_level = 2;
    br_cp.delta = 0.3;
    br_cp.max_individual_dist = 0.6;
    br_cp.bwd_rotations = 4;
    //br_cp.load_from_block(settings_block.get_block("branch_clusterization_params"));

    ClusterizationParams cp;
    cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
    cp.ignore_structure_level = 1;
    cp.delta = 0.3;
    cp.max_individual_dist = ggd.clustering_max_individual_distance;
    cp.bwd_rotations = 4;
    cp.light_importance = 0.8;
    cp.different_types_tolerance = false;
    //cp.load_from_block(settings_block.get_block("tree_clusterization_params"));

    pack_layer(ggd, grove, trees_external, h, packingLayersTrunks, post_voxels,
               tr_cp, 0, 0, true, false, false);

    pack_layer(ggd, grove, trees_external, h, packingLayersBranches, post_voxels,
               br_cp, 1, 1000, true, true, false);

    pack_layer(ggd, grove, trees_external, h, packingLayersTrees, post_voxels,
               cp, 0, 1000, false, false, true);

    transform_all_according_to_root(grove);

    delete (post_voxels);
    delete (post_seeder);
}

void GrovePacker::init()
{
    inited = true;
    BlkManager man;
    man.load_block_from_file("settings.blk",settings_block);   
}

void GrovePacker::pack_grove(GroveGenerationData ggd, GrovePacked &grove, DebugVisualizer &debug,
                             ::Tree *trees_external, Heightmap *h, bool visualize_voxels)
{
    if (!inited)
        init();
    
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
    std::vector<ClusterData> trunks_clusters_base, branches_clusters_base, full_tree_clusters_base;
    std::vector<ClusterData> trunks_clusters, branches_clusters, full_tree_clusters;

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    if (ggd.task & (GenerationTask::MODELS | GenerationTask::SYNTS | GenerationTask::CLUSTERIZE))
    {
        Clusterizer tr_cl;
        tr_cl.set_light(post_voxels);
        tr_cl.get_base_clusters(trees_external, count, 0, trunks_clusters_base);
        if (ggd.task & (GenerationTask::CLUSTERIZE))
        {
            ClusterizationParams tr_cp;
            tr_cp.weights = std::vector<float>{1, 0, 0, 0.0, 0.0};
            tr_cp.ignore_structure_level = 1;
            tr_cp.delta = 0.1;
            tr_cp.light_importance = 0;
            tr_cp.different_types_tolerance = true;
            tr_cp.r_weights = std::vector<float>{0.4, 0, 0, 0.0, 0.0};
            tr_cp.max_individual_dist = 0.4;
            tr_cp.bwd_rotations = 4;
            tr_cl.clusterize(tr_cp, trunks_clusters_base, trunks_clusters);
        }
        else
        {
            trunks_clusters = trunks_clusters_base;
        }
        for (ClusterData &cd : trunks_clusters)
        {
            transform_according_to_root(cd);
        }
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    if (ggd.task & (GenerationTask::MODELS | GenerationTask::SYNTS | GenerationTask::CLUSTERIZE | GenerationTask::BILLBOARDS))
    {

        Clusterizer cl;
        cl.set_light(post_voxels);
        cl.get_base_clusters(trees_external, count, 1, branches_clusters_base);
        if (ggd.task & (GenerationTask::CLUSTERIZE))
        {
            ClusterizationParams cp;
            cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
            cp.ignore_structure_level = 2;
            cp.delta = 0.3;
            cp.max_individual_dist = ggd.clustering_max_individual_distance;
            cp.bwd_rotations = 4;

            cl.clusterize(cp, branches_clusters_base, branches_clusters);
            cl.clusterize(cp, branches_clusters, branches_clusters_base);
            cl.clusterize(cp, branches_clusters_base, branches_clusters);
        }
        else
        {
            branches_clusters = branches_clusters_base;
        }
        for (int i = 0; i < branches_clusters.size(); i++)
        {
            for (::Branch *br : branches_clusters[i].ACDA.originals)
            {
                br->mark_A = i;
            }
        }
    }

    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    if (ggd.task & (GenerationTask::SYNTS))
    {
        SyntheticTreeGenerator stg = SyntheticTreeGenerator(*post_seeder, trunks_clusters, branches_clusters, curGgd);
        stg.generate(trees_external + count, synts, post_voxels);
    }

    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();

    if (ggd.task & (GenerationTask::CLUSTERIZE | GenerationTask::IMPOSTORS | GenerationTask::IMPOSTOR_FULL_GROVE))
    {
        Clusterizer cl2;
        cl2.set_light(post_voxels);
        cl2.get_base_clusters(trees_external, count + synts, 0, full_tree_clusters_base);
        if (ggd.task & (GenerationTask::CLUSTERIZE))
        {
            ClusterizationParams cp;
            cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
            cp.ignore_structure_level = 1;
            cp.delta = 0.3;
            cp.max_individual_dist = ggd.clustering_max_individual_distance;
            cp.bwd_rotations = 4;
            cp.light_importance = 0.8;
            cp.different_types_tolerance = false;
            cl2.clusterize(cp, full_tree_clusters_base, full_tree_clusters);
        }
        else
        {
            full_tree_clusters = full_tree_clusters_base;
        }
    }

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    for (int i = grove.impostors.size(); i < 2; i++)
        grove.impostors.push_back(ImpostorsData());
    if (ggd.task & (GenerationTask::IMPOSTOR_FULL_GROVE))
    {
        logerr("full grove impostors are temporary unavailable");
    }
    if (ggd.task & (GenerationTask::IMPOSTORS))
    {
        ImpostorBaker *ib2 = new ImpostorBaker(ggd.impostor_quality, 0, full_tree_clusters, curGgd.types, &grove.impostors[1]);
        delete (ib2);
    }

    for (int i = grove.clouds.size(); i < 4; i++)
        grove.clouds.push_back(BillboardCloudData());
    if (ggd.task & (GenerationTask::BILLBOARDS))
    {
        BillboardCloudRaw *cloud1 = new BillboardCloudRaw(ggd.bill_1_quality, 1,
                                                          branches_clusters, curGgd.types, &grove.clouds[2]);
        delete (cloud1);
    }
    if (ggd.task & (GenerationTask::BILLBOARDS))
    {
        BillboardCloudRaw *cloud2 = new BillboardCloudRaw(ggd.bill_2_quality, 2,
                                                          branches_clusters, curGgd.types, &grove.clouds[3]);
        delete (cloud2);
    }

    if (ggd.task & (GenerationTask::MODELS))
    {
        std::vector<BranchStructure> instanced_structures;

        for (ClusterData &cd : trunks_clusters)
        {
            pack_cluster(cd, grove, instanced_structures, 0, 0);
        }
        for (ClusterData &cd : branches_clusters)
        {
            pack_cluster(cd, grove, instanced_structures, 1, 1000);
        }
    }
    delete (post_voxels);
    delete (post_seeder);

    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
    /*std::cerr << "Generation took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "[ms]" << std::endl;
    std::cerr << "Main clusterization took " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "[ms]" << std::endl;
    std::cerr << "Syntetic trees generation took " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << "[ms]" << std::endl;
    std::cerr << "Secondary clusterization took " << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count() << "[ms]" << std::endl;
    std::cerr << "Finishing took " << std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count() << "[ms]" << std::endl;
    */
}