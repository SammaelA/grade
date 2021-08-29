#include "branch_clusterization.h"
#include "billboard_cloud.h"
#include "generated_tree.h"
#include "tinyEngine/utility.h"
#include <set>
#include "tinyEngine/utility.h"
#include "GPU_clusterization.h"
int dist_calls = 0;
using namespace glm;
#define DEBUG 0
#define PI 3.14159265f
std::vector<std::pair<float, float>> optimization_quantiles = {
    {0.01, 0.7455902677},
    {0.05, 0.89963103},
    {0.1, 0.9682669918},
    {0.15, 0.9897940455},
    {0.2, 0.996278552},
    {0.25, 0.9984770934},
    {0.3, 0.9992320979}}; //values got by experiments
float distribution[110];
float distribution2[110];
ClusterizationParams clusterizationParams;
long cur_cluster_id = 1;
glm::vec3 Clusterizer::canonical_bbox = glm::vec3(100,60,60);
LightVoxelsCube *test_voxels_cube = nullptr;
struct JSortData
{
    float dist;
    Joint *j1;
    Joint *j2;
    JSortData(float _dist, Joint *_j1, Joint *_j2)
    {
        dist = _dist;
        j1 = _j1;
        j2 = _j2;
    }
    JSortData()
    {
        dist = -1;
        j1 = nullptr;
        j2 = nullptr;
    }
};
struct compare
{
    bool operator()(const JSortData &j1, const JSortData &j2) const
    {
        return j1.dist < j2.dist;
    }
};
ClusterData::ClusterData()
{
    id = cur_cluster_id;
    cur_cluster_id++;
}
inline float AS_branch_min_dist(float part_min_dist, int num_quntile)
{
    return part_min_dist / (1 + optimization_quantiles[num_quntile].first);
}
inline float AS_branch_min_dist(float part_min_dist, float error)
{
    //while measuring distance between branches, we can rotate them and find minimal distance,
    //but we can estimate it from one measurment with particular angle - with this function
    //we use optimization_quantiles to find min_distance so that real distance after rotations
    //will be less than min distance with probability < error
    int nc = -1;
    error = 1 - error;
    for (int i = 0; i < optimization_quantiles.size(); i++)
    {
        if (optimization_quantiles[i].second > error)
        {
            nc = i;
            break;
        }
    }
    if (nc == -1)
        return 0;
    else
        AS_branch_min_dist(part_min_dist, nc);
}

bool dedicated_bbox(Branch *branch, BBox &bbox)
{
    if (!branch || branch->segments.empty())
        return false;
    vec3 a(0, 0, 0);
    vec3 b(0, 0, 0);
    vec3 c;
    a = normalize(branch->joints.back().pos - branch->joints.front().pos);
    for (Joint &j : branch->joints)
    {
        for (Branch *br : j.childBranches)
        {
            b += br->joints.back().pos - br->joints.front().pos;
        }
    }
    if (length(cross(a, b)) < 0.01)
        b = vec3(0, 1, 0);
    b = normalize(b - dot(a, b) * a);
    c = cross(a, b);

    bbox = BillboardCloudRaw::get_bbox(branch, a, b, c);
    return true;
}
Clusterizer::Answer partial_dist(std::vector<int> &jc, std::vector<int> &jp, std::vector<float> &matches, const std::vector<float> &weights)
{
    float num_m = 0.0;
    float num_p = 0.0;
    float denom = 0.0;

    for (int i = 0; i < matches.size(); i++)
    {
        num_m += weights[i] * (2 * matches[i]);
        num_p += weights[i] * jp[i];
        denom += weights[i] * jc[i];
    }

    if (denom < 0.001)
        return Clusterizer::Answer(true, 0, 0);
    num_m /= denom;
    num_p /= denom;
    return Clusterizer::Answer(num_p > 0.9999, num_p - num_m, 1 - num_m);
}
int pass_all_joints(std::vector<int> &jp, Branch *b)
{
    jp[b->level] += b->joints.size();
    for (Joint &j1 : b->joints)
    {
        for (Branch *br : j1.childBranches)
            pass_all_joints(jp, br);
    }
}
bool Clusterizer::match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc,
                                       std::vector<int> &jp, float min, float max)
{
    int sz1 = j1->childBranches.size();
    int sz2 = j2->childBranches.size();
    if (sz1 == 0 || sz2 == 0)
    {
        for (Branch *b : j1->childBranches)
            pass_all_joints(jp, b);
        for (Branch *b : j2->childBranches)
            pass_all_joints(jp, b);
        return true;
    }
    else if (sz1 == 1 && sz2 == 1)
    {
        return match_joints(j1->childBranches.front(), j2->childBranches.front(), matches, jc, jp, min, max);
    }
    else if (sz1 == 1 || sz2 == 1)
    {
        if (sz2 == 1)
        {
            Joint *tmp = j1;
            j1 = j2;
            j2 = tmp;
        }
    }
}
float r_NMSE(Branch *b1, Branch *b2)
{
    double err = 0;
    double sum = 0;
    auto s1 = b1->segments.begin();
    auto s2 = b2->segments.begin();
    while (s1 != b1->segments.end() && s2 != b2->segments.end())
    {
        float rots = MAX(1,MAX(s1->mults.size(),s2->mults.size()));
        for (int i = 0; i < rots;i++)
        {
            float r1 = s1->rel_r_begin*Branch::get_r_mult(2*PI*i/rots,s1->mults);
            float r2 = s2->rel_r_begin*Branch::get_r_mult(2*PI*i/rots,s2->mults);
            err += SQR(r1 - r2);
            sum += SQR(r1) + SQR(r2);
        }
        s1++;
        s2++;
    }
    return err/sum;
}
bool Clusterizer::match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                               float min, float max)
{
    if ((b1->level >= clusterizationParams.ignore_structure_level) &&
        (b2->level >= clusterizationParams.ignore_structure_level))
    {
        return true;
    }
    jp[b1->level] += b1->joints.size() + b2->joints.size();
    if (b1->joints.size() == 1)
    {
        if (b2->joints.size() == 1)
            matches[b1->level]++;
        return true;
    }
    float av_len = 0.5 * (length(b1->joints.back().pos - b1->joints.front().pos) + length(b2->joints.back().pos - b2->joints.front().pos));
    float cur_delta = clusterizationParams.delta * av_len;
    float cur_dist = 0;

    std::multiset<JSortData, compare> distances;
    std::multiset<Joint *> matched_joints;

    for (Joint &j1 : b1->joints)
    {
        for (Joint &j2 : b2->joints)
        {
            float len = length(j1.pos - j2.pos);
            if (len < cur_delta)
            {
                distances.emplace(JSortData(len, &j1, &j2));
            }
        }
    }

    for (Joint &j1 : b1->joints)
    {
        j1.mark_A = 0;
    }
    for (Joint &j2 : b2->joints)
    {
        j2.mark_A = 0;
    }

    auto it = distances.begin();
    while (it != distances.end())
    {
        if (matched_joints.find(it->j1) != matched_joints.end() ||
            matched_joints.find(it->j2) != matched_joints.end())
        {
            it = distances.erase(it);
        }
        else
        {
            matched_joints.emplace(it->j1);
            matched_joints.emplace(it->j2);
            it->j1->mark_A = -1;
            it->j2->mark_A = -1;
            it++;
        }
    }
    //after it only correct matches remained here
    int matches_count = distances.size();
    matches[b1->level] += matches_count;

    if (b1->level == matches.size() - 1) //this is the last branch level, no need to iterate over child branches
        return true;

    for (Joint &j1 : b1->joints)
    {
        if (j1.mark_A >= 0)
        {
            for (Branch *br : j1.childBranches)
            {
                pass_all_joints(jp, br);
            }
        }
    }
    for (Joint &j1 : b2->joints)
    {
        if (j1.mark_A >= 0)
        {
            for (Branch *br : j1.childBranches)
            {
                pass_all_joints(jp, br);
            }
        }
    }
    if (partial_dist(jc, jp, matches, current_data->weights).from > min)
        return false;
    it = distances.begin();
    while (it != distances.end())
    {
        if (!match_child_branches(it->j1, it->j2, matches, jc, jp, min, max))
        {
            //if child branches matching made early exit, it means that min distance limit
            //is already passed. No need to go further.
            return false;
        }
        it++;
    }
    return true;
}
void Clusterizer::get_light(Branch *b, std::vector<float> &light, glm::mat4 &transform)
{
    for (Joint &j : b->joints)
    {
        glm::vec3 ps = transform*glm::vec4(j.pos,1.0f);
        light[b->level] += current_light->get_occlusion(transform*glm::vec4(j.pos,1.0f));
        for (Branch *br : j.childBranches)
        {
            get_light(br,light,transform);
        }
    }
}
Clusterizer::Answer Clusterizer::light_difference(BranchWithData &bwd1, BranchWithData &bwd2)
{
    if (!bwd1.leavesDensity.empty() && !bwd2.leavesDensity.empty())
    {
        if (bwd1.rot_angle < 0)
            bwd1.rot_angle += 2*PI;
        int index = ((int)floor(bwd1.rot_angle/(2*PI) * clusterizationParams.bwd_rotations)) % clusterizationParams.bwd_rotations;
        float res = bwd1.leavesDensity[index]->NMSE(bwd2.leavesDensity[0]);
        return Answer(true, res, res);
    }
    else
        return Answer(true,0,0);
    if (!current_light)
        return Answer(true,0,0);
    std::vector<float> _l11(current_data->light_weights.size(),0);
    std::vector<float> _l12(current_data->light_weights.size(),0);
    std::vector<float> _l21(current_data->light_weights.size(),0);
    std::vector<float> _l22(current_data->light_weights.size(),0);
    glm::mat4 t1 = bwd1.transform; 
    glm::mat4 t2 = bwd2.transform; 
    get_light(bwd1.b,_l11,t1);//real b1
    get_light(bwd2.b,_l22,t2);//real b2
    get_light(bwd1.b,_l12,t2);//b1 placed instead of b2
    get_light(bwd2.b,_l21,t1);//b2 placed instead of b1
    double l11 = 0, l12 = 0, l21 = 0 ,l22 = 0;
    for (int i=0;i<current_data->light_weights.size();i++)
    {
        l11 += current_data->light_weights[i]*_l11[i];
        l12 += current_data->light_weights[i]*_l12[i];
        l21 += current_data->light_weights[i]*_l21[i];
        l22 += current_data->light_weights[i]*_l22[i];
    }
    double res = (abs(l12 - l11) + abs(l21 - l22))/(l11 + l12 + l21 + l22);
    return Answer(true, res, res);
}
Clusterizer::Answer Clusterizer::dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    if ((b1->type_id != b2->type_id && !clusterizationParams.different_types_tolerance)|| 
         b1->level != b2->level)
        return Answer(true,1000,1000);
    std::vector<int> joint_counts(bwd1.joint_counts);
    std::vector<int> joint_passed(joint_counts.size(), 0);
    std::vector<float> matches(joint_counts.size());
    float light_importance = clusterizationParams.light_importance;
    for (int i = 0; i < joint_counts.size(); i++)
    {
        joint_counts[i] += bwd2.joint_counts[i];
    }
    bool exact = match_joints(b1, b2, matches, joint_counts, joint_passed, min/(1 - light_importance), max);
    Answer part_answer = partial_dist(joint_counts, joint_passed, matches, current_data->weights);

    if (exact)
    {
        Answer light_answer = light_difference(bwd1, bwd2);
        Answer res = light_answer*(light_importance) + part_answer*(1 - light_importance);
        float r_importance = clusterizationParams.r_weights[CLAMP(b1->level,0,clusterizationParams.r_weights.size()-1)];
        if (r_importance > 0)
        {
            float r_diff = r_importance*r_NMSE(b1,b2);
            res.from += r_diff;
            res.to += r_diff;
        }
        return res;
    }
    else
        return part_answer;
}
Clusterizer::Answer Clusterizer::dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
}
Clusterizer::Answer Clusterizer::dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
{
    bwd1.rot_angle = 0;
    Answer fast_answer = dist_simple(bwd1, bwd2, min, max);
    if (!fast_answer.exact)
    {
        return fast_answer;
    }
    min = min > fast_answer.from ? fast_answer.from : min;
    int N = 5;
    int iterations = 4;
    float min_dist = fast_answer.from;
    float min_phi = 0;
    float base_step = 2 * PI / N;
    vec3 axis = bwd1.b->joints.back().pos - bwd1.b->joints.front().pos;
    logerr("axis = %f %f %f",axis.x,axis.y,axis.z);
    mat4 rot = rotate(mat4(1.0f), base_step, axis);
    for (int i = 1; i <= N; i++)
    {
        bwd1.b->transform(rot);
        bwd1.rot_angle = i * base_step;
        float md = dist_simple(bwd1, bwd2, min, max).from;
        if (md < min_dist)
        {
            min_dist = md;
            min_phi = i * base_step;
        }
    }
    mat4 rot_to_base = rotate(mat4(1.0f), min_phi, axis);
    for (int i = 0; i < iterations; i++)
    {
        base_step = 0.5 * base_step;

        rot = rot_to_base * rotate(mat4(1.0f), base_step, axis);
        bwd1.b->transform(rot);
        bwd1.rot_angle = min_phi + base_step;
        float md_plus = dist_simple(bwd1, bwd2, min, max).from;

        rot = rotate(mat4(1.0f), -2.0f * base_step, axis);
        bwd1.b->transform(rot);
        bwd1.rot_angle = min_phi - base_step;
        float md_minus = dist_simple(bwd1, bwd2, min, max).from;
        if (md_plus < md_minus && md_plus < min_dist)
        {
            min_dist = md_plus;
            min_phi += base_step;
            rot_to_base = rotate(mat4(1.0f), 2.0f * base_step, axis);
        }
        else if (md_minus < md_plus && md_minus < min_dist)
        {
            min_dist = md_minus;
            min_phi -= base_step;
            rot_to_base = rotate(mat4(1.0f), 0.0f, axis);
        }
        else
        {
            rot_to_base = rotate(mat4(1.0f), 1.0f * base_step, axis);
        }
    }
    rot_to_base = rot_to_base * rotate(mat4(1.0f), -min_phi, axis);
    bwd1.b->transform(rot_to_base);
    if (data)
        data->rotation = min_phi;
    return Answer(true, min_dist, min_dist);
}
Clusterizer::Answer Clusterizer::dist_trunc(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max,
                                            DistData *data)
{
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    if (data)
        data->rotation = 0;
    if (b1->joints.size() != b2->joints.size())
        return Answer(true,1000,1000);
    if (b1->joints.size()<2)
        return Answer(true,0,0);
    auto it1 = b1->joints.begin();
    auto it2 = b2->joints.begin();
    glm::vec3 prev1 = it1->pos;
    glm::vec3 prev2 = it2->pos;
    it1++;
    it2++;
    float av_len = 0.5 * (length(b1->joints.back().pos - b1->joints.front().pos) + length(b2->joints.back().pos - b2->joints.front().pos));
    float cur_delta = clusterizationParams.delta * av_len;
    
    while (it1!=b1->joints.end() && it2!=b2->joints.end())
    {
        if (length(it1->pos-it2->pos) > cur_delta)
            return Answer(true,1000,1000);
        it1++;
        it2++;
    }
}
Clusterizer::Answer Clusterizer::dist(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
{
    dist_calls++;
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    if ((b1->type_id != b2->type_id && !clusterizationParams.different_types_tolerance) || b1->level != b2->level)
        return Answer(true,1000,1000);
    return dist_Nsection(bwd1, bwd2, min, max, data);
}
void Clusterizer::get_base_clusters(Tree &t, int layer, std::vector<ClusterData> &base_clusters)
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
                base_clusters.back().IDA.type_ids.push_back(b.type_id);
                base_clusters.back().IDA.tree_ids.push_back(t.id);
                base_clusters.back().IDA.centers_par.push_back(b.center_par);
                base_clusters.back().IDA.centers_self.push_back(b.center_self);
                base_clusters.back().IDA.transforms.push_back(glm::mat4(1.0f));
                base_clusters.back().ACDA.originals.push_back(nullptr);
                //since we delete full trees right after packing the in clusters originals now mean nothing
                base_clusters.back().ACDA.ids.push_back(b.self_id);
                base_clusters.back().ACDA.rotations.push_back(0);
        }
    }
}
void Clusterizer::calc_joints_count(Branch *b, std::vector<int> &counts)
{
    for (Joint &j : b->joints)
    {
        counts[b->level]++;
        for (Branch *br : j.childBranches)
        {
            calc_joints_count(br, counts);
        }
    }
}
void Clusterizer::get_base_clusters(Tree *t, int count, int layer, std::vector<ClusterData> &base_clusters)
{
    ClusterizationTmpData data;
    current_data = &data;
    Cluster::currentClusterizer = this;
    
    for (int i = 0; i < count; i++)
    {
        int prev_n = base_clusters.size();
        get_base_clusters(t[i], layer, base_clusters);
        debugl(3, " added %d branches from tree %d\n", base_clusters.size() - prev_n, i);
    }
    current_data = nullptr;
}
void Clusterizer::set_branches(std::vector<ClusterData> &base_clusters)
{
    int i = 0;
    for (ClusterData &cd : base_clusters)
    {
        if (!cd.base || (cd.ACDA.originals.empty()) || cd.IDA.transforms.empty())
            continue;
            BBox bbox;
            Branch &b = *(cd.base);
            Branch *nb = current_data->branchHeap.new_branch();
            nb->deep_copy(&b, current_data->branchHeap, &current_data->leafHeap);
            if (dedicated_bbox(nb, bbox))
            {
                glm::vec3 cbb = canonical_bbox;
                mat4 rot_inv(vec4(bbox.a, 0), vec4(bbox.b, 0), vec4(bbox.c, 0), vec4(0, 0, 0, 1));
                mat4 rot = inverse(rot_inv);
                vec3 sc_vert = vec3(MAX((1/cbb.x) * bbox.sizes.x,MAX( (1/cbb.y) * bbox.sizes.y, (1/cbb.z) * bbox.sizes.z)));
                float r_transform = 1/sc_vert.x;
                mat4 SC = scale(mat4(1.0f), sc_vert);
                mat4 SC_inv = inverse(SC);
                vec3 base_joint_pos = vec4(b.joints.front().pos, 1.0f);
                mat4 transl = translate(mat4(1.0f), -1.0f * base_joint_pos);
                rot = SC_inv * rot * transl;
                nb->transform(rot, r_transform);
                current_data->branches.push_back(BranchWithData(&b, nb, i, MAX_BRANCH_LEVELS, current_data->branches.size(), inverse(rot)));
            }
            i++;
    }
}
void Clusterizer::set_light(LightVoxelsCube *_light)
{
    current_light = _light;
}
ClusterData Clusterizer::extract_data(std::vector<ClusterData> &base_clusters, Clusterizer::Cluster &cl)
{
    ClusterData cd;
    cd.base = cl.prepare_to_replace(base_clusters, cd.IDA, cd.ACDA, cd.id);
    return cd;
}
Clusterizer::~Clusterizer()
{
    if (current_data)
        delete current_data;
}
void Clusterizer::clusterize(ClusterizationParams &params, std::vector<ClusterData> &base_clusters, 
                             std::vector<ClusterData> &clusters, bool need_only_ddt)
{
    ClusterizationTmpData *data = new ClusterizationTmpData();
    data->light_weights = params.light_weights;
    data->weights = params.weights;

    current_data = data;
    Cluster::currentClusterizer = this;
    clusterizationParams = params;

    set_branches(base_clusters);
    if (clusterizationParams.different_types_tolerance == false && !test_voxels_cube)
    {
        test_voxels_cube = new LightVoxelsCube(current_data->branches.front().leavesDensity.front());
    }
    dist_calls = 0;
    prepare_ddt();
    if (!need_only_ddt)
    {
        current_data->Ddg.make_base_clusters(current_data->branches);
        current_data->Ddg.make(20, clusterizationParams.min_clusters);
        for (int c_num : current_data->Ddg.current_clusters)
        {
            clusters.push_back(extract_data(base_clusters, current_data->Ddg.clusters[c_num]));
        }
    }
    for (auto &b : current_data->branches)
    {
        b.clear();
    }
    current_data->branches.clear();
}
void Clusterizer::visualize_clusters(DebugVisualizer &debug, bool need_debug)
{
    std::vector<ClusterData> base_clusters;
    std::vector<ClusterData> _clusters;
    ClusterizationParams params;
    clusterize(params, base_clusters, _clusters);

    if (!need_debug)
        return;
    std::vector<Branch *> branches;
    int k = 0;
    for (int S : current_data->Ddg.current_clusters)
    {
        current_data->Ddg.clusters[S].to_branch_data(branches);
        for (int i = 0; i < branches.size(); i++)
        {
            debug.add_branch_debug(branches[i], vec3(1, 1, 1), vec3(k, -100, 0), -1);
            debug.add_branch_debug(branches[i], vec3(1, 1, 1), vec3(k, 0, 0), 2);
            debug.add_branch_debug(branches[i], vec3(1, 1, 1), vec3(k, 100, 0), 3);
        }
        branches.clear();
        k += 100;
    }
}
Clusterizer::ClusterDendrogramm::Dist
Clusterizer::ClusterDendrogramm::get_P_delta(int n, std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta)
{
    debugl(1, "get P delta\n");
    Dist md(-1, -1, 1000);
    int k = n > current_clusters.size() ? current_clusters.size() : n;
    k = n > current_clusters.size() / 3 ? n : current_clusters.size() / 3;
    int n1 = 0;
    auto i = current_clusters.begin();
    auto j = current_clusters.rbegin();
    while (n1 < k)
    {
        if (*i == *j)
        {
        }
        else
        {
            float distance = clusters[*i].ward_dist(&(clusters[*j]), md.d);
            if (distance < md.d && distance > 0.001)
            {
                md.d = distance;
                md.U = *i;
                md.V = *j;
            }
        }
        i++;
        j++;
        n1++;
    }
    delta = md.d;
    debugl(1, "P delta delta %f\n", delta);
    for (int u : current_clusters)
    {
        for (int v : current_clusters)
        {
            if (u != v)
            {
                int k = -1,l = -1;
                if (clusters[u].branch)
                    k = clusters[u].branch->b->type_id;
                if (clusters[v].branch)
                    l = clusters[v].branch->b->type_id;
                float distance = clusters[u].ward_dist(&(clusters[v]), delta);
                if (distance <= delta)
                {
                    P_delta.push_back(Dist(u, v, distance));
                    if (distance < md.d)
                    {
                        md.d = distance;
                        md.U = u;
                        md.V = v;
                    }
            }
            }
        }
    }
    debugl(1, "P_delta size = %d md = (%d %d %f)\n", P_delta.size(), md.U, md.V, md.d);
    return md;
}
void Clusterizer::ClusterDendrogramm::make(int n, int clusters_num)
{
    std::list<Dist> P_delta;
    float delta;
    Dist min = get_P_delta(n, current_clusters, P_delta, delta);
    for (int i = 1; i < size; i++)
    {
        if (P_delta.empty())
            min = get_P_delta(n, current_clusters, P_delta, delta);
        else if (min.U < 0)
        {
            for (Dist &d : P_delta)
            {
                if (d.U == d.V)
                    debugl(0, "error in P_delta %d\n", d.U);
                if (d.d < min.d)
                {
                    min.U = d.U;
                    min.V = d.V;
                    min.d = d.d;
                }
            }
        }
        if (min.d > 1000*clusterizationParams.max_individual_dist || current_clusters.size() <= clusters_num)
        {
            //logerr("breaking clusterization %f %d %d %d",min.d, (int)clusterizationParams.max_individual_dist, (int)(current_clusters.size()), clusters_num);
            break;
            //makes no sense to merge clusters with maximum distance between them.
        }
        current_clusters.remove(min.U);
        current_clusters.remove(min.V);
        clusters.push_back(Cluster(&(clusters[min.U]), &(clusters[min.V])));
        int W = clusters.size() - 1;
        auto dit = P_delta.begin();
        while (dit != P_delta.end())
        {
            Dist &d = *dit;
            if (d.U == min.U || d.V == min.U || d.U == min.V || d.V == min.V)
                dit = P_delta.erase(dit);
            else
                dit++;
        }
        for (int S : current_clusters)
        {
            float d = clusters[W].ward_dist(&(clusters[S]), delta);
            if (d < delta)
            {
                P_delta.push_back(Dist(W, S, d));
            }
        }
        current_clusters.push_back(W);
        min = Dist(-1, -1, 1000);
        int sum = 0;
        for (int S : current_clusters)
        {
            sum += clusters[S].size;
        }
    }
    int sum = 0;
    for (int S : current_clusters)
    {
        debugl(1, "cluster %d size = %d\n", S, clusters[S].size);
        sum += clusters[S].size;
    }
    float comp = (float)size/current_clusters.size();
    debugl(17, "%d %d %f clusters elements compression\n", current_clusters.size(), size, comp);
}

Clusterizer::Answer Clusterizer::get_dist(BranchWithData &bwd1, BranchWithData &bwd2, DistData *data)
{
    auto p = current_data->ddt.get(bwd1.id,bwd2.id);
    if (data)
        *data = p.second;
    return p.first;
}
void Clusterizer::prepare_ddt()
{
    current_data->ddt.create(current_data->branches.size());
    for (int i = 0; i < current_data->branches.size(); i++)
    {
        current_data->pos_in_table_by_id.emplace(current_data->branches[i].original->self_id,i);
    }
    GPUClusterizationHelper gpuch;
    gpuch.prepare_ddt(current_data->branches,current_data->ddt,clusterizationParams);
    return;
    for (int i = 0; i < current_data->branches.size(); i++)
    {
        for (int j = 0; j < current_data->branches.size(); j++)
        {
            Answer a;
            DistData d;
            if (i == j)
            {
                a.exact = true;
                a.from = 0;
                a.to = 0;
                d.dist = 0;
                d.rotation = 0;
            }
            else if (j < i)
            {
                auto p = current_data->ddt.get(j,i);
                a = p.first;
                d = p.second;
            }
            else
            {
                a = dist(current_data->branches[i],current_data->branches[j],clusterizationParams.max_individual_dist,0,&d);
            }
            current_data->ddt.set(i,j,a,d);
        }
    }
}

void ClusterizationParams::load_from_block(Block *b)
{
    if (!b)
        return;
    
    bwd_rotations = b->get_int("bwd_rotations",bwd_rotations);
    delta = b->get_double("delta",delta);
    light_importance = b->get_double("light_importance",light_importance);
    voxels_size_mult = b->get_double("voxels_size_mult",voxels_size_mult);
    ignore_structure_level = b->get_int("ignore_structure_level",ignore_structure_level);
    min_clusters = b->get_int("min_clusters",min_clusters);
    max_individual_dist = b->get_double("max_individual_dist",max_individual_dist);
    different_types_tolerance = b->get_bool("different_types_tolerance",different_types_tolerance);
    b->get_arr("weights",weights, true);
    b->get_arr("light_weights",light_weights, true);
    b->get_arr("r_weights",r_weights, true);
}

Clusterizer::BranchWithData::~BranchWithData()
{
}

void voxelize_branch(Branch *b, LightVoxelsCube *light, int level_to);
Clusterizer::BranchWithData::BranchWithData(Branch *_original, Branch *_b, int _base_cluster_id, int levels, int _id, glm::mat4 _transform)
{
    base_cluster_id = _base_cluster_id;
    original = _original;
    b = _b;
    id = _id;
    transform = _transform;
    for (int i = 0; i < levels; i++)
        joint_counts.push_back(0);
    if (b)
        calc_joints_count(b, joint_counts);

    glm::vec3 axis = b->joints.back().pos - b->joints.front().pos;
    glm::mat4 rot = glm::rotate(glm::mat4(1.0f), 2 * PI / clusterizationParams.bwd_rotations, axis);

    for (int i = 0; i < clusterizationParams.bwd_rotations; i++)
    {
        b->transform(rot);
        leavesDensity.push_back(new LightVoxelsCube(
            glm::vec3(0.5f*canonical_bbox.x,0,0),
            glm::vec3(0.5f*canonical_bbox.x,canonical_bbox.y, canonical_bbox.z),
            1 / clusterizationParams.voxels_size_mult, 1));
        //set_occlusion(b, leavesDensity.back());
        voxelize_branch(b, leavesDensity.back(), 1);
    }
}

void Clusterizer::BranchWithData::clear()
{
    for (int i = 0; i < leavesDensity.size(); i++)
    {
        if (leavesDensity[i])
        {
            delete leavesDensity[i];
            leavesDensity[i] = nullptr;
        }
    }
    leavesDensity.clear();
}

void Clusterizer::BranchWithData::set_occlusion(Branch *b, LightVoxelsCube *light)
{
    for (Joint &j : b->joints)
    {
        if (j.leaf)
            light->set_occluder_trilinear(j.pos, 1);
        for (Branch *br : j.childBranches)
        {
            set_occlusion(br, light);
        }
    }
}

void voxelize_branch(Branch *b, LightVoxelsCube *light, int level_to)
{
    glm::vec3 second_vec = glm::vec3(b->plane_coef.x,b->plane_coef.y,b->plane_coef.z);
    for (Segment &s : b->segments)
    {
        glm::vec3 dir = s.end - s.begin;
        glm::vec3 n = glm::normalize(glm::cross(dir, second_vec));
        float len = length(dir);
        float v = (1/3.0)*PI*(SQR(s.rel_r_begin) + s.rel_r_begin*s.rel_r_end + SQR(s.rel_r_end));
        float R = 0.5*(s.rel_r_begin + s.rel_r_end);
        int samples = MIN(5 + 0.25*v, 50);

        //approximate partial cone with cylinder
        //uniformly distributed points in cylinder
        for (int i = 0; i<samples; i++)
        {
            float phi = urand(0, 2*PI);
            float h = urand(0, 1);
            float r = R*sqrtf(urand(0, 1));
            glm::vec3 nr = glm::rotate(glm::mat4(1),phi,dir)*glm::vec4(n,0);
            glm::vec3 sample = s.begin + h*dir + r*nr;
            light->set_occluder_trilinear(sample, 25/samples);
        } 
    }
    if (b->level < level_to)
    {
        for (Joint &j : b->joints)
        {
            for (Branch *br : j.childBranches)
            {
                voxelize_branch(br, light, level_to);
            }
        }
    }
}
void Clusterizer::get_current_ddt(Clusterizer::DistDataTable &ddt)
{
    if (current_data)
    {
        ddt.copy(current_data->ddt);
    }
}

std::map<int,int> Clusterizer::get_id_pos_map()
{
    if (current_data)
    {
        return current_data->pos_in_table_by_id;
    }
    else
    {
        return std::map<int,int>();
    }
}

void Clusterizer::set_default_clustering_params(ClusterizationParams &params, ClusteringStep step)
{
    if (step == ClusteringStep::TRUNKS)
    {
        ClusterizationParams tr_cp;
        tr_cp.weights = std::vector<float>{1, 0, 0, 0.0, 0.0};
        tr_cp.ignore_structure_level = 1;
        tr_cp.delta = 0.1;
        tr_cp.light_importance = 0;
        tr_cp.different_types_tolerance = true;
        tr_cp.r_weights = std::vector<float>{0.4, 0, 0, 0.0, 0.0};
        tr_cp.max_individual_dist = 0.6;
        tr_cp.bwd_rotations = 4;
        params = tr_cp;
    }
    else if (step == ClusteringStep::BRANCHES)
    {
        ClusterizationParams br_cp;
        br_cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
        br_cp.ignore_structure_level = 2;
        br_cp.delta = 0.3;
        br_cp.max_individual_dist = 0.6;
        br_cp.bwd_rotations = 4;
        params = br_cp;
    }
    else if (step == ClusteringStep::TRUNKS)
    {
        ClusterizationParams cp;
        cp.weights = std::vector<float>{5000, 800, 40, 0.0, 0.0};
        cp.ignore_structure_level = 1;
        cp.delta = 0.3;
        cp.max_individual_dist = 0.7;
        cp.bwd_rotations = 4;
        cp.light_importance = 0.8;
        cp.different_types_tolerance = false;
        params = cp;
    }
}