#include "clustering.h"
#include "helpers.h"
#include "dummies.h"
#include "hierarcial_clustering.h"
#include "structural_similarity.h"
#include "hasing.h"
#include "impostor_metric.h"
#include "clustering.h"
#include "visualize_clusters.h"
#include "pyclustering.h"
#include "GPU_impostor_metric.h"

ClusteringStrtegy cStrategy = ClusteringStrtegy::Merge;

using namespace glm;
int cur_cluster_id = 0;
ClusterData::ClusterData()
{
    id = cur_cluster_id;
    cur_cluster_id++;
}
void Clusterizer2::prepare(Block &settings)
{
    Block &def_b = get_default_block();
    std::string c_helper_name = def_b.get_string("clustering_helper","impostor");
    c_helper_name = settings.get_string("clustering_helper",c_helper_name);

    std::string c_base_name = def_b.get_string("clustering_base","hierarcial");
    c_base_name = settings.get_string("clustering_base",c_base_name);

    if (c_helper_name == "impostor")
        clusteringHelper = new ImpostorClusteringHelper();
    else if (c_helper_name == "gpu_impostor")
        clusteringHelper = new GPUImpostorClusteringHelper();
    else if (c_helper_name == "structural_similarity_cpu")
        clusteringHelper = new CPUSSClusteringHelper();
    else if (c_helper_name == "structural_similarity_gpu")
        clusteringHelper = new GPUSSClusteringHelper();
    else if (c_helper_name == "hash_ddt")
        clusteringHelper = new DDTHashBasedClusteringHelper();
    else if (c_helper_name == "hash_simple")
        clusteringHelper = new SimpleHashBasedClusteringHelper();
    else
    {
        logerr("given unknown clustering helper name %s",c_helper_name);
        clusteringHelper = new ImpostorClusteringHelper();
    }
    if (c_base_name == "hierarcial")
        clusteringBase = new HierarcialClusteringBase();
    else if (c_base_name == "py_kmeans")
        clusteringBase = new KmeansPyClusteringBase();
    else if (c_base_name == "py_xmeans")
        clusteringBase = new XmeansPyClusteringBase();
    else if (c_base_name == "dbscan")
        clusteringBase = new XmeansPyClusteringBase();
    else if (c_base_name == "optics")
        clusteringBase = new XmeansPyClusteringBase();
    else
    {
        logerr("given unknown clustering base name %s",c_base_name);
        clusteringBase = new HierarcialClusteringBase();
    }
}
Clusterizer2::~Clusterizer2()
{
    if (clusteringHelper)
        delete clusteringHelper;
    if (clusteringBase)
        delete clusteringBase;
}
void Clusterizer2::get_base_clusters(Block &settings, Tree *t, int count, int layer, std::vector<ClusterData> &base_clusters,
                                     ClusteringContext *ctx)
{
    for (int i = 0; i < count; i++)
    {
        int prev_n = base_clusters.size();
        get_base_clusters(settings, t[i], layer, base_clusters, ctx);
        debugl(3, " added %d branches from tree %d\n", base_clusters.size() - prev_n, i);
    }
}
BranchClusteringData *Clusterizer2::convert_branch(Block &settings, Branch *base, ClusteringContext *ctx)
{
    BaseBranchClusteringData bbcd;
    BBox bbox;
    Branch &b = *(base);
    if (get_dedicated_bbox(base, bbox))
    {
        mat4 transform = mat4(1.0f);
        glm::vec3 cbb = canonical_bbox();
        mat4 rot_inv(vec4(bbox.a, 0), vec4(bbox.b, 0), vec4(bbox.c, 0), vec4(0, 0, 0, 1));
        mat4 rot = inverse(rot_inv);
        vec3 sc_vert = vec3(MAX((1 / cbb.x) * bbox.sizes.x, MAX((1 / cbb.y) * bbox.sizes.y, (1 / cbb.z) * bbox.sizes.z)));
        float r_transform = 1 / sc_vert.x;
        mat4 SC = scale(mat4(1.0f), sc_vert);
        mat4 SC_inv = inverse(SC);
        vec3 base_joint_pos = vec4(b.joints.front().pos, 1.0f);
        mat4 transl = translate(mat4(1.0f), -1.0f * base_joint_pos);
        rot = SC_inv * rot * transl;
        transform = inverse(rot);
        bbcd.r_transform = r_transform;
        bbcd.transform = transform;
        bbcd.can_be_center = true;
    }
    BranchClusteringData *br = clusteringHelper->convert_branch(settings, base, ctx, bbcd);
    br->set_base(bbcd);
    return br;
}
void Clusterizer2::get_base_clusters(Block &settings, Tree &t, int layer, std::vector<ClusterData> &base_clusters,
                                     ClusteringContext *ctx)
{
    if (!t.valid)
        return;
    if (layer < 0 || layer >= t.branchHeaps.size() || t.branchHeaps[layer]->branches.size() == 0)
    {
        return;
    }
    else
    {
        for (Branch &b : t.branchHeaps[layer]->branches)
        {
            base_clusters.push_back(ClusterData());
            base_clusters.back().base = &b;
            base_clusters.back().base_pos = 0;
            base_clusters.back().IDA.type_ids.push_back(b.type_id);
            base_clusters.back().IDA.tree_ids.push_back(t.id);
            base_clusters.back().IDA.centers_par.push_back(b.center_par);
            base_clusters.back().IDA.centers_self.push_back(b.center_self);
            base_clusters.back().IDA.transforms.push_back(glm::mat4(1.0f));
            base_clusters.back().ACDA.originals.push_back(nullptr);
            //since we delete full trees right after packing the in clusters originals now mean nothing
            base_clusters.back().ACDA.ids.push_back(b.self_id);
            base_clusters.back().ACDA.rotations.push_back(0);
            base_clusters.back().ACDA.clustering_data.push_back(convert_branch(settings, &b, ctx));
        }
    }
}
void Clusterizer2::prepare_branches(Block &settings, std::vector<ClusterData> &base_clusters,
                                    std::vector<BranchClusteringData *> &branches)
{
    /*for (auto &c : base_clusters)
    {
        for (BranchClusteringData *b : c.ACDA.clustering_data)
        {
            if (b)
                branches.push_back(b);
        }
    }*/
    if (cStrategy == ClusteringStrtegy::Merge)
    {
        for (auto &c : base_clusters)
        {
            branches.push_back(c.ACDA.clustering_data[c.base_pos]);
            branches.back()->base_cluster_id = c.id;
            logerr("cluster id %d %d %d",c.id, c.base_pos, branches.back());
        }
    }

    for (int i = 0; i < base_clusters.size(); i++)
    {
        tmpData.pos_in_table_by_id.emplace(base_clusters[i].id, i);
        tmpData.pos_in_table_by_branch_id.emplace(base_clusters[i].base->self_id, i);
    }
}
void Clusterizer2::clusterize(Block &settings, std::vector<ClusterData> &base_clusters, std::vector<ClusterData> &clusters,
                              ClusteringContext *ctx, bool need_save_full_data)
{
    tmpData = ClusterizationTmpData();
    std::vector<BranchClusteringData *> branches;
    prepare_branches(settings, base_clusters, branches);

    IntermediateClusteringData *ICD = clusteringHelper->prepare_intermediate_data(settings, branches, ctx);
    ICD->branches = branches;
    ICD->elements_count = branches.size();
    std::vector<ClusteringBase::ClusterStruct> cluster_result;
    clusteringBase->clusterize(settings, ICD, cluster_result);

    prepare_result(settings, base_clusters, clusters, branches, ctx, cluster_result);
    //if (current_clustering_step == ClusteringStep::BRANCHES)
    //    visualize_clusters(settings, branches, cluster_result, ctx, "clusters",128,128);
    if (need_save_full_data)
    {
        fcd = new FullClusteringData();
        fcd->base_clusters = &base_clusters;
        fcd->clusters = cluster_result;
        fcd->result_clusters = &clusters;
        fcd->id = ICD;
        fcd->ctx = ctx;
        fcd->settings = &settings;
        fcd->pos_in_table_by_id = tmpData.pos_in_table_by_id;
        fcd->pos_in_table_by_branch_id = tmpData.pos_in_table_by_branch_id;
    }
    else
    {
        delete ICD;
        if (cStrategy == ClusteringStrtegy::Merge)
        {
            //we do not need clustering data for non-central branches
            for (auto &cl : clusters)
            {
                for (int i = 0; i < cl.ACDA.clustering_data.size(); i++)
                {
                    if (i == cl.base_pos)
                        continue;
                    if (cl.ACDA.clustering_data[i])
                    {
                        clusteringHelper->clear_branch_data(cl.ACDA.clustering_data[i], ctx);
                        cl.ACDA.clustering_data[i] = nullptr;
                    }
                }
            }
        }
    }

    tmpData = ClusterizationTmpData();
}
void Clusterizer2::prepare_result(Block &settings, std::vector<ClusterData> &base_clusters, std::vector<ClusterData> &result_clusters,
                                  std::vector<BranchClusteringData *> &branches, ClusteringContext *ctx,
                                  std::vector<ClusteringBase::ClusterStruct> &result)
{
    if (cStrategy == ClusteringStrtegy::Merge)
    {
        for (auto *b : branches)
        {
            logerr("center %d", b->base_cluster_id);
        }
        for (auto &str : result)
        {
            result_clusters.emplace_back();
            auto &res_cluster = result_clusters.back();
            auto &IDA = res_cluster.IDA;
            auto &ACDA = res_cluster.ACDA;
            BranchClusteringData *center = branches[str.center];
            logerr("center %d from %d id %d",str.center, branches.size(),center->base_cluster_id);
            auto it = tmpData.pos_in_table_by_id.find(center->base_cluster_id);
            if (it == tmpData.pos_in_table_by_id.end())
            {
                logerr("cannot find central base cluster by it's id %d", center->base_cluster_id);
                return;
            }
            res_cluster.base = base_clusters[it->second].base;
            res_cluster.id = base_clusters[it->second].id;
            glm::mat4 base_transform_inv = glm::inverse(center->transform);
            for (auto &p : str.members)
            {
                BranchClusteringData *base_bcd = branches[p.first];
                glm::mat4 rot = glm::rotate(glm::mat4(1.0f),p.second.rot,glm::vec3(1,0,0));
                logerr("base bcd first %d %d %d",base_bcd, center, p.first);
                glm::mat4 tr = (base_bcd->transform) * rot * base_transform_inv;

                it = tmpData.pos_in_table_by_id.find(base_bcd->base_cluster_id);
                if (it == tmpData.pos_in_table_by_id.end())
                {
                    logerr("cannot find central base cluster by it's id %d", base_bcd->base_cluster_id);
                    return;
                }
                ClusterData &base = base_clusters[it->second];
                if (base_bcd == center)
                {
                    res_cluster.base_pos = ACDA.originals.size() + base.base_pos;
                }
                if (!base.is_valid())
                {
                    logerr("failed to merge clusters. Base cluster is not valid and will be ignored");
                }
                else
                {
                    int sz = base.ACDA.originals.size();
                    for (int i = 0; i < sz; i++)
                    {
                        IDA.transforms.push_back(base.IDA.transforms[i] * tr);
                        IDA.centers_par.push_back(base.IDA.centers_par[i]);
                        IDA.centers_self.push_back(base.IDA.centers_self[i]);
                        IDA.type_ids.push_back(base.IDA.type_ids[i]);
                        IDA.tree_ids.push_back(base.IDA.tree_ids[i]);
                        ACDA.rotations.push_back(base.ACDA.rotations[i]);
                        ACDA.originals.push_back(base.ACDA.originals[i]);
                        ACDA.ids.push_back(base.ACDA.ids[i]);
                        ACDA.clustering_data.push_back(base.ACDA.clustering_data[i]);
                    }
                }
            }
        }
    }
}