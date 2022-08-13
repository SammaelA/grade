#include "structural_similarity.h"
#include "GPU_clusterization.h"
#include <set>

LightVoxelsCube *current_light = nullptr;
ClassicStructureSimilarityParams clusterizationParams;
namespace structsim
{
    Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data);
}
BranchClusteringData *StructuralSimilarityClusteringHelper::convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, 
                                                                           BaseBranchClusteringData &data)
{
    clusterizationParams.load(&settings);
    BranchWithData *bwd = new BranchWithData(clusterizationParams, base, MAX_BRANCH_LEVELS, 0, data.transform, data.r_transform);
    
    return bwd;
}
void StructuralSimilarityClusteringHelper::clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx)
{
    if (base)
    {
        base->clear();
        delete base;
    }
}

IntermediateClusteringData *CPUSSClusteringHelper::prepare_intermediate_data(Block &settings, 
                                                                             std::vector<BranchClusteringData *> branches,
                                                                             ClusteringContext *ctx)
{
    clusterizationParams.load(&settings);
    current_light = nullptr;//TODO: to resurrect this algorithm you need to calculate current_light right here 

    IntermediateClusteringDataDDT *data = new IntermediateClusteringDataDDT();

    data->branches = branches;
    data->ddt.create(branches.size());
    data->elements_count = branches.size();
    std::vector<BranchWithData *> real_branches;
    for (int i = 0; i < branches.size(); i++)
    {
        real_branches.push_back(dynamic_cast<BranchWithData *>(branches[i]));
        if (!real_branches.back())
        {
            logerr("CPUSSClusteringHelper error: wrong type of BranchClusteringData");
            return data;
        }
    }
    for (int i = 0; i < branches.size(); i++)
    {
        for (int j = 0; j < branches.size(); j++)
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
                auto p = data->ddt.get(j,i);
                a = p.first;
                d = p.second;
            }
            else
            {
                a = structsim::dist(*(real_branches[i]),*(real_branches[j]),clusterizationParams.max_individual_dist,0,&d);
                logerr("dist[%d %d] = %f",i,j,a.from);
            }
            data->ddt.set(i,j,a,d);
        }
    }

    return data;
}
IntermediateClusteringData *GPUSSClusteringHelper::prepare_intermediate_data(Block &settings, 
                                                                             std::vector<BranchClusteringData *> branches,
                                                                             ClusteringContext *ctx)
{
    clusterizationParams.load(&settings);
    current_light = nullptr;//TODO: to resurrect this algorithm you need to calculate current_light right here 

    IntermediateClusteringDataDDT *data = new IntermediateClusteringDataDDT();

    data->branches = branches;
    data->ddt.create(branches.size());
    data->elements_count = branches.size();
    std::vector<BranchWithData *> real_branches;
    for (int i = 0; i < branches.size(); i++)
    {
        real_branches.push_back(dynamic_cast<BranchWithData *>(branches[i]));
        if (!real_branches.back())
        {
            logerr("CPUSSClusteringHelper error: wrong type of BranchClusteringData");
            return data;
        }
    }
    GPUClusterizationHelper gpuch;
    gpuch.prepare_ddt(real_branches,data->ddt,clusterizationParams);
}
namespace structsim
{
    
    using namespace glm;
    int dist_calls;
    std::vector<std::pair<float, float>> optimization_quantiles = {
    {0.01, 0.7455902677},
    {0.05, 0.89963103},
    {0.1, 0.9682669918},
    {0.15, 0.9897940455},
    {0.2, 0.996278552},
    {0.25, 0.9984770934},
    {0.3, 0.9992320979}}; //values got by experiments

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
bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                               float min, float max);
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
Answer partial_dist(std::vector<int> &jc, std::vector<int> &jp, std::vector<float> &matches, const std::vector<float> &weights)
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
        return Answer(true, 0, 0);
    num_m /= denom;
    num_p /= denom;
    return Answer(num_p > 0.9999, num_p - num_m, 1 - num_m);
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
bool match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc,
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
bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
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
    if (partial_dist(jc, jp, matches, clusterizationParams.weights).from > min)
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
void get_light(Branch *b, std::vector<float> &light, glm::mat4 &transform)
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
Answer light_difference(BranchWithData &bwd1, BranchWithData &bwd2)
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
    std::vector<float> _l11(clusterizationParams.light_weights.size(),0);
    std::vector<float> _l12(clusterizationParams.light_weights.size(),0);
    std::vector<float> _l21(clusterizationParams.light_weights.size(),0);
    std::vector<float> _l22(clusterizationParams.light_weights.size(),0);
    glm::mat4 t1 = bwd1.transform; 
    glm::mat4 t2 = bwd2.transform; 
    get_light(bwd1.b,_l11,t1);//real b1
    get_light(bwd2.b,_l22,t2);//real b2
    get_light(bwd1.b,_l12,t2);//b1 placed instead of b2
    get_light(bwd2.b,_l21,t1);//b2 placed instead of b1
    double l11 = 0, l12 = 0, l21 = 0 ,l22 = 0;
    for (int i=0;i<clusterizationParams.light_weights.size();i++)
    {
        l11 += clusterizationParams.light_weights[i]*_l11[i];
        l12 += clusterizationParams.light_weights[i]*_l12[i];
        l21 += clusterizationParams.light_weights[i]*_l21[i];
        l22 += clusterizationParams.light_weights[i]*_l22[i];
    }
    double res = (abs(l12 - l11) + abs(l21 - l22))/(l11 + l12 + l21 + l22);
    return Answer(true, res, res);
}
Answer dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
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
    Answer part_answer = partial_dist(joint_counts, joint_passed, matches, clusterizationParams.weights);

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
Answer dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
}
Answer dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
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
Answer dist_trunc(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max,
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
Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
{
    dist_calls++;
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    if ((b1->type_id != b2->type_id && !clusterizationParams.different_types_tolerance) || b1->level != b2->level)
        return Answer(true,1000,1000);
    return dist_Nsection(bwd1, bwd2, min, max, data);
}
}
