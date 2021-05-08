#include "synthetic_trees_generator.h"
#include <iostream>
SyntheticTreeGenerator::SyntheticTreeGenerator(Seeder &_seeder, std::vector<ClusterData> &_trunks_clusters,
                                               std::vector<ClusterData> &_branches_clusters, GroveGenerationData &_ggd) : seeder(_seeder),
                                                                                                                          ggd(_ggd),
                                                                                                                          trunks_clusters(_trunks_clusters),
                                                                                                                          branches_clusters(_branches_clusters)
{
}
void SyntheticTreeGenerator::get_mean_and_stddev(std::vector<float> &values, double &mean, double &stddev)
{
    if (values.empty())
    {
        logerr("stat error - empty sample given");
        return;
    }
    if (values.size() == 1)
    {
        debugl(5, "stat warning - cannot calculate stddev by 1 sample");
        mean = values.front();
        stddev = 0;
        return;
    }
    mean = 0;
    stddev = 0;
    for (float &val : values)
    {
        mean += val;
    }
    mean /= values.size();
    for (float &val : values)
    {
        stddev += SQR(val - mean);
    }
    stddev = sqrt(stddev / (values.size() - 1));
}
BranchStat SyntheticTreeGenerator::get_branch_stat(ClusterData &cd)
{
    BranchStat bs;
    bs.valid = true;
    std::vector<float> phi_s, psi_s, r_s, rots;
    std::map<int, int> type_id_and_count;

    for (int i = 0; i < cd.ACDA.originals.size(); i++)
    {
        Branch *b = cd.ACDA.originals[i];
        if (type_id_and_count.empty())
        {
            type_id_and_count.emplace(b->type_id, 1);
        }
        else
        {
            auto it = type_id_and_count.find(b->type_id);
            if (it == type_id_and_count.end())
            {
                type_id_and_count.emplace(b->type_id, 1);
            }
            else
            {
                it->second++;
            }
        }
        rots.push_back(cd.ACDA.rotations[i]);
        glm::vec3 dir = b->joints.size() > 2 ? b->joints.back().pos - b->joints.front().pos : glm::vec3(1, 1, 1);
        r_s.push_back(length(dir));
        dir = glm::normalize(dir);
        psi_s.push_back(acos(dir.y));
        if (abs(dir.x) < 1e-9)
            phi_s.push_back(dir.z > 0 ? PI / 2 : -PI / 2);
        else
            phi_s.push_back(atan(dir.z / dir.x));
    }

    std::vector<double> values;
    std::vector<double> weights;
    for (auto &pr : type_id_and_count)
    {
        values.push_back(pr.first);
        weights.push_back(pr.second);
    }
    bs.typeStat = new DiscreteGeneral(values, weights);
    double mean, dev;
    float sample_size = phi_s.size() + 1;
    float ssq = pow(sample_size, 1.5);
    //calculate distribution of branch parameters. Increase deviation if sample size is small
    get_mean_and_stddev(phi_s, mean, dev);
    bs.transformStat.phi = new Normal(mean, dev + 2 * PI / (ssq * ssq));
    get_mean_and_stddev(psi_s, mean, dev);
    bs.transformStat.psi = new Normal(mean, 0.5 * dev + PI / (3 * ssq * ssq));
    get_mean_and_stddev(r_s, mean, dev);
    bs.transformStat.r = new Normal(mean, dev + 0.25 * mean / ssq);
    get_mean_and_stddev(rots, mean, dev);
    bs.transformStat.rot_angle = new Normal(mean, dev + PI / (2 * ssq));
    return bs;
}
DiscreteGeneral *SyntheticTreeGenerator::get_child_branches_cluster_stat(ClusterData &cd)
{
    std::map<int, int> cluster_id_and_count;
    int child_branches_count = 0;
    for (int i = 0; i < cd.ACDA.originals.size(); i++)
    {
        Branch *b = cd.ACDA.originals[i];
        for (Joint &j : b->joints)
        {
            for (Branch *chb : j.childBranches)
            {
                child_branches_count++;
                auto it = cluster_id_and_count.find(chb->mark_A);
                if (it == cluster_id_and_count.end())
                {
                    cluster_id_and_count.emplace(chb->mark_A, 1);
                }
                else
                {
                    it->second += 1;
                }
            }
        }
    }
    std::vector<double> values;
    std::vector<double> weights;
    for (auto &pr : cluster_id_and_count)
    {
        values.push_back(pr.first);
        weights.push_back(pr.second);
    }
    if (values.empty())
        return nullptr;
    else
        return new DiscreteGeneral(values, weights);
}
void SyntheticTreeGenerator::get_existance_stat(ClusterData &cd, std::vector<DiscreteGeneral *> &branchExistanceStat)
{
    const int max_child_branches = 4;
    std::vector<std::vector<double>> cbs;
    for (int i = 0; i < cd.ACDA.originals.size(); i++)
    {
        Branch *b = cd.ACDA.originals[i];
        int k = 0;
        for (Joint &j : b->joints)
        {
            if (k == cbs.size())
            {
                cbs.push_back(std::vector<double>(max_child_branches, 0));
            }
            cbs[k][MIN(max_child_branches - 1, j.childBranches.size())]++;
            k++;
        }
    }
    std::vector<double> values;
    for (int i = 0; i < max_child_branches; i++)
    {
        values.push_back(i);
    }
    for (auto &weights : cbs)
    {
        branchExistanceStat.push_back(new DiscreteGeneral(values, weights));
    }
}
void SyntheticTreeGenerator::add_shadow(Branch *b, glm::mat4 &transform)
{
    if (!b)
        return;
    for (Joint &j : b->joints)
    {
        glm::vec3 ps = transform*glm::vec4(j.pos,1);
        voxels->set_occluder(ps, powf(4 - b->level, 2));
        for (Branch *br : j.childBranches)
        {
            add_shadow(br, transform);
        }
    }
}
void SyntheticTreeGenerator::add_planar_shadow(Tree &t)
{
    seeder.add_tree_shadow(t);
}
void SyntheticTreeGenerator::generate(Tree *_trees, int count, LightVoxelsCube *_voxels)
{
    voxels = _voxels;
    collect_statistics();
    for (int i = 0; i < count; i++)
    {
        SyntheticTree t;
        make_synt_tree(t);
        synt_tree_to_real(t, _trees[i]);
        add_planar_shadow(_trees[i]);
    }
}
void SyntheticTreeGenerator::collect_statistics()
{

    int i = 0;
    std::vector<double> weights;
    std::vector<double> types;
    for (ClusterData &cd : trunks_clusters)
    {
        stat.rootStats.push_back(RootStat());
        stat.rootStats.back().selfBranchStat = get_branch_stat(cd);
        stat.rootStats.back().childBranchesClusterStat = get_child_branches_cluster_stat(cd);
        get_existance_stat(cd, stat.rootStats.back().branchExistanceStat);
        if (!stat.rootStats.back().childBranchesClusterStat)
        {
            std::vector<double> v = {0};
            std::vector<double> w = {0};
            stat.rootStats.back().childBranchesClusterStat = new DiscreteGeneral(v, w);
            weights.push_back(0);
            types.push_back(i);
            stat.rootStats.back().child_branches_count = 0;
        }
        else
        {
            weights.push_back(cd.IDA.transforms.size());
            types.push_back(i);
            stat.rootStats.back().child_branches_count = cd.base->joints.size();
        }
        i++;
    }
    stat.rootClusterStat = new DiscreteGeneral(types, weights);

    for (ClusterData &cd : branches_clusters)
    {
        stat.branchStats.push_back(get_branch_stat(cd));
    }
}
void SyntheticTreeGenerator::synt_tree_to_real(SyntheticTree &synt, Tree &t)
{
    for (int i = 0; i < synt.trunk_instance_data.type_ids.size(); i++)
    {
        synt.trunk->IDA.type_ids.push_back(synt.trunk_instance_data.type_ids[i]);
        synt.trunk->IDA.centers_par.push_back(synt.trunk_instance_data.centers_par[i]);
        synt.trunk->IDA.centers_self.push_back(synt.trunk_instance_data.centers_self[i]);
        synt.trunk->IDA.transforms.push_back(synt.trunk_instance_data.transforms[i]);
    }
    for (int j = 0; j < synt.branches.size(); j++)
    {
        for (int i = 0; i < synt.branches_instance_data[j].type_ids.size(); i++)
        {
            synt.branches[j]->IDA.type_ids.push_back(synt.branches_instance_data[j].type_ids[i]);
            synt.branches[j]->IDA.centers_par.push_back(synt.branches_instance_data[j].centers_par[i]);
            synt.branches[j]->IDA.centers_self.push_back(synt.branches_instance_data[j].centers_self[i]);
            synt.branches[j]->IDA.transforms.push_back(synt.branches_instance_data[j].transforms[i]);
        }
    }

    t.leaves = new LeafHeap();
    BranchHeap *bh = new BranchHeap();
    t.branchHeaps.push_back(bh);
    BranchHeap *bh1 = new BranchHeap();
    t.branchHeaps.push_back(bh1);
    t.root = bh->new_branch();
    t.root->norecursive_copy(synt.trunk->base, *bh1, t.leaves);
    t.root->transform(synt.trunk_instance_data.transforms.front());
    t.root->type_id = synt.trunk_instance_data.type_ids.front();
    t.root->center_par = synt.trunk_instance_data.centers_par.front();
    t.root->center_self = synt.trunk_instance_data.centers_self.front();
    t.type = &(ggd.types[t.root->type_id]);
    int k = 0;
    for (Joint &j : t.root->joints)
    {
        for (int i = 0; i < synt.branches.size(); i++)
        {
            if (synt.joints_ns[i] == k)
            {
                Branch *ch = bh1->new_branch();
                ch->deep_copy(synt.branches[i]->base, *bh1, t.leaves);
                ch->transform(synt.branches_instance_data[i].transforms.front());
                j.childBranches.push_back(ch);
            }
        }
        k++;
    }
}
glm::mat4 SyntheticTreeGenerator::get_transform(Branch *base, glm::vec3 pos, BranchStat &stat)
{
    glm::vec3 dir = base->joints.back().pos - base->joints.front().pos;
    if (glm::length(dir) < 1e-4)
        return glm::mat4(0);
    float scale = stat.transformStat.r->get();
    scale /= glm::length(dir);
    dir = glm::normalize(dir);
    float psi_s = stat.transformStat.psi->get() - acos(dir.y);
    float phi_s = stat.transformStat.phi->get();
    if (abs(dir.x) < 1e-9)
        phi_s -= (dir.z > 0 ? PI / 2 : -PI / 2);
    else
        phi_s -= atan(dir.z / dir.x);
    if (base->level == 0)
    {
        //TODO replace with something better
        phi_s *= 0.2;
        psi_s *= 0.2;
    }
    float self_rot = stat.transformStat.rot_angle->get();
    glm::mat4 transf = glm::translate(
        glm::rotate(
            glm::rotate(
                glm::rotate(
                    glm::scale(
                        glm::translate(glm::mat4(1.0f), pos),
                        glm::vec3(scale)),
                    self_rot, dir),
                phi_s, glm::vec3(0, 1, 0)),
            psi_s, glm::vec3(1, 0, 0)),
        -(base->joints.front().pos));
    glm::vec4 ps = glm::vec4(base->joints.front().pos, 1);
    glm::vec4 p = transf * ps;
    if (p.x != p.x || p.y != p.y || p.z != p.z)
    {
        logerr("matrix is NAN %f %f %f) pos = (%f %f %f) scale = %f, self_rot = %f, dir = (%f %f %f),phi_s = %f %f",
               ps.x, ps.y, ps.z, pos.x, pos.y, pos.z, scale, self_rot, dir.x, dir.y, dir.z, psi_s, phi_s);
    }
    return transf;
}
float calc_occ(Branch *b, glm::mat4 &transform, LightVoxelsCube *voxels, int max_level)
{
    if (!b || !voxels || b->level > max_level)
        return 0;
    if (b->level == max_level)
    {
        glm::vec3 ps = transform * glm::vec4(b->joints.front().pos, 1);
        return voxels->get_occlusion(ps);
    }
    float occ = 0;
    for (Joint &j : b->joints)
    {
        glm::vec3 ps = transform * glm::vec4(j.pos, 1);
        occ += voxels->get_occlusion(ps);
        for (Branch *br : j.childBranches)
        {
            occ += calc_occ(br,transform,voxels,max_level);
        }
    }
    return occ;
}
float calc_occ_non_rec(Branch *b, glm::mat4 &transform, LightVoxelsCube *voxels)
{
    if (!b || !voxels)
        return 0;
    float occ = 0;
    for (Joint &j : b->joints)
    {
        glm::vec3 ps = transform * glm::vec4(j.pos, 1);
        occ += voxels->get_occlusion(ps);
    }
    return occ;
}
void SyntheticTreeGenerator::make_synt_tree(SyntheticTree &synt)
{
    int max_general_iters = 1;
    int max_tries = pow(2,ggd.synts_precision);
    int max_root_iters = pow(2,ggd.synts_precision);
    float threshold = 0.1;
    float max_occlusion = 1e6;
    float eps = 0.01;
    struct Bdata
    {
        int br_cl = -1;
        glm::mat4 transform = glm::mat4(0);
        float occlusion = 1e9;
    };
    struct RootBdata
    {
        int br_cl = -1;
        glm::mat4 transform = glm::mat4(0);
        float occlusion = 1e9;
        glm::vec3 pos;
    };
    RootBdata bd;

    for (int k = 0; k < max_general_iters; k++)
    {
        std::vector<RootBdata> min_root_bdata;
        float min_root_occ = max_occlusion;
        synt = SyntheticTree();
        for (int i = 0; i < max_root_iters; i++)
        {
            std::vector<Seed> seeds;
            seeder  .choose_places_for_seeds(1, seeds);
            if (seeds.empty())
                continue;

            glm::vec3 pos = glm::vec3(seeds.front().pos.x, 0, seeds.front().pos.y);
            pos.y = seeder.heightmap->get_height(pos);
            int root_cl;
            RootStat *root_stat = nullptr;
            root_cl = CLAMP(stat.rootClusterStat->get(), 0, trunks_clusters.size() - 1);
            root_stat = &(stat.rootStats[root_cl]);

            if (root_stat->child_branches_count <= 0 || trunks_clusters[root_cl].base->joints.size() == 0)
                continue;
            
            RootBdata bdata;
            bdata.br_cl = root_cl;
            bdata.pos = pos;
            bdata.transform = get_transform(trunks_clusters[root_cl].base, pos, root_stat->selfBranchStat);
            bdata.occlusion = calc_occ_non_rec(trunks_clusters[root_cl].base,bdata.transform,voxels);
            if (bdata.occlusion < (1 - eps) * min_root_occ)
            {
                min_root_occ = bdata.occlusion;
                min_root_bdata.clear();
                min_root_bdata.push_back(bdata);
            }
            else if (bdata.occlusion < (1 + eps) * min_root_occ)
            {
                min_root_occ = MIN(min_root_occ, bdata.occlusion);
                min_root_bdata.push_back(bdata);
            }
        }

        if (!min_root_bdata.empty() && min_root_occ < max_occlusion)
        {
            int bdind = urandi(0,min_root_bdata.size());
            bd = min_root_bdata[bdind];
            RootStat &root_stat = stat.rootStats[bd.br_cl];
            synt.trunk = &(trunks_clusters[bd.br_cl]);
            synt.trunk_instance_data.type_ids.push_back(root_stat.selfBranchStat.typeStat->get());
            synt.trunk_instance_data.centers_par.push_back(ggd.pos);
            synt.trunk_instance_data.centers_self.push_back(bd.pos);
            synt.trunk_instance_data.transforms.push_back(bd.transform);
        }
        else
        {
            logerr("synthetic tree creation failed - impossible to find source root with child branches");
            continue;
        }

        int i = 0;
        int synt_child_branches_count = 0;
        int root_cl = bd.br_cl;
        RootStat *root_stat = &(stat.rootStats[bd.br_cl]);
        glm::vec3 root_pos = bd.pos;
        for (Joint &base_j : trunks_clusters[root_cl].base->joints)
        {
            int chb_count = root_stat->branchExistanceStat[MIN(root_stat->branchExistanceStat.size(), i)]->get();
            if (i == trunks_clusters[root_cl].base->joints.size() - 1)
                chb_count = MAX(1, chb_count);
            glm::vec3 bpos = synt.trunk_instance_data.transforms.front() * glm::vec4(base_j.pos, 1);
            for (int j = 0; j < chb_count; j++)
            {
                int child_branches_added = 0;

                std::vector<Bdata> min_bdata;
                float min_occl = max_occlusion;
                for (int k = 0; k < max_tries; k++)
                {
                    Bdata bdata;
                    int c = root_stat->childBranchesClusterStat->get();
                    int br_cl = CLAMP(c, 0, branches_clusters.size() - 1);
                    Branch *base = branches_clusters[br_cl].base;
                    auto b_stat = stat.branchStats[br_cl];

                    glm::mat4 transform = get_transform(base, bpos, b_stat);
                    float occ = calc_occ(base,transform,voxels,2);

                    bdata.br_cl = br_cl;
                    bdata.transform = transform;
                    bdata.occlusion = occ;
                    if (bdata.occlusion < (1 - eps) * min_occl)
                    {
                        min_occl = bdata.occlusion;
                        min_bdata.clear();
                        min_bdata.push_back(bdata);
                    }
                    else if (bdata.occlusion < (1 + eps) * min_occl)
                    {
                        min_occl = MIN(min_occl, bdata.occlusion);
                        min_bdata.push_back(bdata);
                    }
                }
                if (min_occl < max_occlusion)
                {
                    int rnd = urandi(0, min_bdata.size());
                    Bdata rnd_min = min_bdata[rnd];
                    add_shadow(branches_clusters[rnd_min.br_cl].base, rnd_min.transform);
                    if (child_branches_added == 0)
                        synt.joints_ns.push_back(i);
                    int c = root_stat->childBranchesClusterStat->get();
                    int br_cl = rnd_min.br_cl;
                    auto b_stat = stat.branchStats[br_cl];
                    synt.branches.push_back(&(branches_clusters[br_cl]));
                    synt.branches_instance_data.push_back(InstanceDataArrays());
                    synt.branches_instance_data.back().centers_par.push_back(root_pos);
                    synt.branches_instance_data.back().centers_self.push_back(bpos);
                    int id = synt.trunk_instance_data.type_ids.back();
                    synt.branches_instance_data.back().type_ids.push_back(id);
                    synt.branches_instance_data.back().transforms.push_back(rnd_min.transform);
                    child_branches_added++;
                    synt_child_branches_count++;
                }
            }
            i++;
        }
        if (synt_child_branches_count > threshold * trunks_clusters[root_cl].base->joints.size())
        {
            return;
        }
    }
}
FullStat::~FullStat()
{
#define DEL(a) \
    if (a)     \
        delete (a);
    DEL(rootClusterStat);
    for (auto &st : rootStats)
    {
        DEL(st.childBranchesClusterStat);
        for (auto *st2 : st.branchExistanceStat)
        {
            DEL(st2);
        }
        DEL(st.selfBranchStat.typeStat);
        DEL(st.selfBranchStat.transformStat.phi);
        DEL(st.selfBranchStat.transformStat.psi);
        DEL(st.selfBranchStat.transformStat.r);
        DEL(st.selfBranchStat.transformStat.rot_angle);
    }
    for (auto &st : branchStats)
    {
        DEL(st.typeStat);
        DEL(st.transformStat.phi);
        DEL(st.transformStat.psi);
        DEL(st.transformStat.r);
        DEL(st.transformStat.rot_angle);
    }
}