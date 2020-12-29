#include "branch_clusterization.h"
#include "billboard_cloud.h"
#include "generated_tree.h"
#include <set>
using namespace glm;
#define DEBUG 0
#define PI 3.14159265f
std::vector<std::pair<float, float>> optimization_quantiles = {
    {0.01,0.7455902677},
    {0.05,0.89963103},
    {0.1,0.9682669918},
    {0.15,0.9897940455},
    {0.2,0.996278552},
    {0.25,0.9984770934},
    {0.3,0.9992320979}
};//values got by experiments
int distribution[110];
std::vector<float> Clusterizer::weights = std::vector<float>{0.6, 0.25, 0.1, 0.035, 0.015};
float Clusterizer::delta = 0.25;
struct JSortData
{
    float dist;
    Joint *j1;
    Joint *j2;
    JSortData(float _dist, Joint *_j1,Joint *_j2)
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
inline float AS_branch_min_dist(float part_min_dist, int num_quntile)
{
    return part_min_dist/(1 + optimization_quantiles[num_quntile].first);
}
inline float AS_branch_min_dist(float part_min_dist, float error)
{
    //while measuring distance between branches, we can rotate them and find minimal distance, 
    //but we can estimate it from one measurment with particular angle - with this function
    //we use optimization_quantiles to find min_distance so that real distance after rotations
    //will be less than min distance with probability < error
    int nc = -1;
    error = 1 - error;
    for (int i=0;i<optimization_quantiles.size();i++)
    {
        if (optimization_quantiles[i].second > error)
        {
            nc = i;
            break;
        }
    }
    if (nc == -1)
        return 0;
    else AS_branch_min_dist(part_min_dist,nc);
}

void transform_branch(Branch *b, mat4 transform)
{
    b->base_seg_n = b->joints.size();
    //in clusterization process we use this variable
    //for joints count is branch and child branches to speed-up 
    //calculations 
    for (Segment &s : b->segments)
    {
        s.begin = transform*vec4(s.begin,1.0f);
        s.end = transform*vec4(s.end,1.0f);
    }
    for (Joint &j : b->joints)
    {

        j.pos = transform*vec4(j.pos,1.0f);
        for (Branch *br : j.childBranches)
        {
            transform_branch(br,transform);
            b->base_seg_n += br->base_seg_n;
        }
    }
}
bool dedicated_bbox(Branch *branch, BillboardCloud::BBox &bbox)
{
    if (!branch || false && branch->dead || branch->segments.empty())
        return false;
    vec3 a(0,0,0);
    vec3 b(0,0,0);
    vec3 c;
    a = normalize(branch->joints.back().pos - branch->joints.front().pos);
    for (Joint &j : branch->joints)
    {
        for (Branch *br : j.childBranches)
        {
            b += br->joints.back().pos - br->joints.front().pos;
        }
    }
    if (length(cross(a,b))<0.01)
        b = vec3(0,1,0);
    b = normalize(b - dot(a,b)*a);
    c = cross(a,b);

    bbox = BillboardCloud::get_bbox(branch,a,b,c);
    return true;
}
Clusterizer::Answer partial_dist(std::vector<int> &jc, std::vector<int> &jp, std::vector<float> &matches,  const std::vector<float> &weights)
{
    float num_m=0.0;
    float num_p=0.0;
    float denom = 0.0;
    for (int i=0;i<matches.size();i++)
    {
        num_m += weights[i]*(2*matches[i]);
        num_p += weights[i]*jp[i];
        denom +=weights[i]*jc[i];
    }
    if (denom < 0.001)
        return Clusterizer::Answer(true,0,0);
    num_m /= denom;
    num_p /= denom;
    //fprintf(stderr,"numm nump %f %f\n",num_m, num_p);
    return Clusterizer::Answer(num_p > 0.9999,num_p - num_m, 1 - num_m);
}
int pass_all_joints(std::vector<int> &jp, Branch *b)
{
    jp[b->level] += b->joints.size();
    for (Joint &j1 : b->joints)
    {
        for (Branch *br : j1.childBranches)
            pass_all_joints(jp,br);
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
            pass_all_joints(jp,b);
        for (Branch *b : j2->childBranches)
            pass_all_joints(jp,b);
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
bool Clusterizer::match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                               float min, float max)
{
    jp[b1->level] += b1->joints.size() + b2->joints.size();
    if (b1->joints.size() == 1)
    {
        if (b2->joints.size() == 1)
            matches[b1->level]++;
        return true;
    }
    float av_len = 0.5*(length(b1->joints.back().pos - b1->joints.front().pos) + length(b2->joints.back().pos - b2->joints.front().pos));
    float cur_delta = delta*av_len;
    float cur_dist = 0;

    std::multiset<JSortData,compare> distances;
    std::multiset<Joint *> matched_joints;

    for (Joint &j1 : b1->joints)
    {
        for (Joint &j2 : b2->joints)
        {
            float len = length(j1.pos - j2.pos);
            if (len<cur_delta)
            {
                distances.emplace(JSortData(len, &j1, &j2));
            }
        }
    }

    for (Joint &j1 : b1->joints)
    {
        j1.max_branching = 0;
    }
    for (Joint &j2 : b2->joints)
    {
        j2.max_branching = 0;
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
            it->j1->max_branching = -1;
            it->j2->max_branching = -1;
            it++;
        }
    }
    //after it only correct matches remained here
    int matches_count = distances.size();
    matches[b1->level] += matches_count;

    if (b1->level == matches.size() - 1)//this is the last branch level, no need to iterate over child branches
        return true;
    
    for (Joint &j1 : b1->joints)
    {
        if (j1.max_branching>=0)
        {
            for (Branch *br : j1.childBranches)
            {
                pass_all_joints(jp,br);
            }
        }
    }
    for (Joint &j1 : b2->joints)
    {
        if (j1.max_branching>=0)
        {   
            for (Branch *br : j1.childBranches)
            {
                pass_all_joints(jp,br);
            }
        }
    }
    if (partial_dist(jc,jp,matches,weights).from > min)
        return false;
    it = distances.begin();
    while (it != distances.end())
    {
        if (!match_child_branches(it->j1,it->j2,matches,jc,jp,min,max))
        {
            //if child branches matching made early exit, it means that min distance limit
            //is already passed. No need to go further.
            return false;
        }
        it++;
    }
    return true;

}
Clusterizer::Answer Clusterizer::dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    Answer part_answer;
    std::vector<int> joint_counts(bwd1.joint_counts);
    std::vector<int> joint_passed(joint_counts.size(),0);
    std::vector<float> matches(joint_counts.size());
    for (int i = 0; i < joint_counts.size(); i++)
    {
        joint_counts[i] += bwd2.joint_counts[i];
    }
    bool exact = match_joints(b1, b2, matches, joint_counts, joint_passed, min, max);
    part_answer = partial_dist(joint_counts, joint_passed, matches, weights);
    return part_answer;
}
Clusterizer::Answer Clusterizer::dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
    /*float min_dist = 2;
    int min_i = 0;
    float dr = 2*PI/360;
    vec3 axis = bwd1.b->joints.back().pos - bwd1.b->joints.front().pos;
    mat4 rot = rotate(mat4(1.0f),dr,axis);
    //fprintf(stderr,"axis = %f %f %f\n",axis.x,axis.y,axis.z);
    float dists[360];
    for (int i=1;i<=360;i++)
    {

        Branch *b1 = bwd1.b;
        b1->transform(rot);
        Branch *b2 = bwd2.b;
        Answer part_answer;
        std::vector<int> joint_counts(bwd1.joint_counts);
        std::vector<float> matches(joint_counts.size());
        for (int i=0;i<joint_counts.size();i++)
        {
            joint_counts[i] += bwd2.joint_counts[i];
        }
        part_answer.exact = match_joints(b1,b2,matches,joint_counts,min,max);
        part_answer.from = partial_dist(joint_counts,matches,weights);
        part_answer.to = part_answer.exact ? part_answer.from : 1;
        if (part_answer.from<min_dist)
        {
            min_dist = part_answer.from;
            min_i = i;
        }
        dists[i % 360] = part_answer.from;
    }
    int loc_mins_count = 0;
    for (int i=0;i<360;i++)
    {
        if ((dists[i] < dists[(i+1)%360]) && (dists[i] < dists[(i-1)%360]))
            loc_mins_count++;
    }
    if (loc_mins_count > 35)
        loc_mins_count = 35;
    rot = rotate(mat4(1.0f),min_i*dr,axis);
    Branch *b1 = bwd1.b;
        b1->transform(rot);
        Branch *b2 = bwd2.b;
        Answer part_answer;
        std::vector<int> joint_counts(bwd1.joint_counts);
        std::vector<float> matches(joint_counts.size());
        for (int i=0;i<joint_counts.size();i++)
        {
            joint_counts[i] += bwd2.joint_counts[i];
        }
        part_answer.exact = match_joints(b1,b2,matches,joint_counts,min,max);
        part_answer.from = partial_dist(joint_counts,matches,weights);
        part_answer.to = part_answer.exact ? part_answer.from : 1;
    //fprintf(stderr,"final dist = %f\n",part_answer.from);
    return part_answer;*/
}
Clusterizer::Answer Clusterizer::dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
{
    Answer fast_answer = dist_simple(bwd1,bwd2,min,max);
    if (!fast_answer.exact)
    {
        return fast_answer;
    }
    min = min > fast_answer.from ? fast_answer.from : min;
    int N = 5;
    int iterations = 4;
    float min_dist = fast_answer.from;
    float min_phi = 0;
    float base_step = 2*PI/N;
    vec3 axis = bwd1.b->joints.back().pos - bwd1.b->joints.front().pos;
    mat4 rot = rotate(mat4(1.0f),base_step,axis);
    for (int i=1;i<=N;i++)
    {
        bwd1.b->transform(rot);
        float md = dist_simple(bwd1,bwd2,min,max).from;
        if (md<min_dist)
        {
            min_dist = md;
            min_phi = i*base_step;
        }
    }
    mat4 rot_to_base = rotate(mat4(1.0f),min_phi,axis);
    for (int i=0;i<iterations;i++)
    {
        base_step = 0.5*base_step;

        rot = rot_to_base*rotate(mat4(1.0f),base_step,axis);
        bwd1.b->transform(rot);
        float md_plus = dist_simple(bwd1,bwd2,min,max).from;

        rot = rotate(mat4(1.0f),-2.0f*base_step,axis);
        bwd1.b->transform(rot);
        float md_minus = dist_simple(bwd1,bwd2,min,max).from;
        if (md_plus < md_minus && md_plus < min_dist)
        {
            min_dist = md_plus;
            min_phi += base_step;
            rot_to_base = rotate(mat4(1.0f),2.0f*base_step,axis);
        }
        else if (md_minus < md_plus && md_minus < min_dist)
        {
            min_dist = md_minus;
            min_phi -= base_step;
            rot_to_base = rotate(mat4(1.0f),0.0f,axis);
        }
        else
        {
            rot_to_base = rotate(mat4(1.0f),1.0f*base_step,axis);
        }
    }
    rot_to_base = rot_to_base*rotate(mat4(1.0f),-min_phi,axis);
    bwd1.b->transform(rot_to_base);
    if (data)
        data->rotation = min_phi;
    return Answer(true,min_dist,min_dist);
}
Clusterizer::Answer Clusterizer::dist(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max, DistData *data)
{
    return dist_Nsection(bwd1,bwd2,min,max,data);
}
bool Clusterizer::set_branches(Tree &t, int layer)
{
    if (layer<0 || layer>t.branchHeaps.size() || t.branchHeaps[layer]->branches.size() == 0)
    {
        return false;
    }
    else
    {
        int i = 0;
        for (Branch &b : t.branchHeaps[layer]->branches)
        {
            BillboardCloud::BBox bbox;
            Branch *nb = branchHeap.new_branch();
            nb->deep_copy(&b,branchHeap);
            if (dedicated_bbox(nb,bbox)) 
            {
            
                mat4 rot_inv(vec4(bbox.a,0),vec4(bbox.b,0),vec4(bbox.c,0),vec4(0,0,0,1));
                mat4 rot = inverse(rot_inv);
                vec3 base_joint_pos = rot*vec4(b.joints.front().pos,1.0f);
                mat4 transl = translate(mat4(1.0f),-1.0f*base_joint_pos);
                mat4 SC = scale(mat4(1.0f), vec3(0.01f*bbox.sizes.x,0.05f*bbox.sizes.y,0.05f*bbox.sizes.z));
                mat4 SC_inv = inverse(SC);
                rot = SC_inv*transl*rot;
                transform_branch(nb,rot);
                branches.push_back(BranchWithData(&b,nb,t.params.max_depth(),i,inverse(rot)));
                i++;
                
            }
        }    
    }   
}
void Clusterizer::calc_joints_count(Branch *b, std::vector<int> &counts)
{
    if (b->dead)
        return;
    for (Joint &j : b->joints)
    {
        counts[b->level]++;
        for (Branch *br : j.childBranches)
        {
            calc_joints_count(br,counts);
        }
    }
}
bool Clusterizer::set_branches(Tree *t, int count, int layer)
{
    for (int i=0;i<count;i++)
    {
        int prev_n = branches.size();
        set_branches(t[i],layer);
        fprintf(stderr," added %d branches from tree %d\n", branches.size() - prev_n, i);
    }
}
void Clusterizer::visualize_clusters(DebugVisualizer &debug, bool need_debug)
{
    ClusterDendrogramm Ddg;
    Ddg.make_base_clusters(branches);
    Ddg.make(20);
    fprintf(stderr,"dist decreasion distribution:");
    for (int i=0;i<100;i++)
    {
        fprintf(stderr," %d,",distribution[i]);
    }
    fprintf(stderr,"\n");
    if (!need_debug)
        return;
    std::vector<Branch *> branches;
    int k = 0;
    for (int S : Ddg.current_clusters)
    {
        Ddg.clusters[S].to_branch_data(branches);
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
Clusterizer::ClusterDendrogramm::get_P_delta(int n,std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta)
{
    fprintf(stderr, "get P delta\n");
    Dist md(-1,-1,1000);
    int k = n > current_clusters.size() ? current_clusters.size() : n;
    k = n > current_clusters.size()/3 ? n : current_clusters.size()/3;
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
            float distance = clusters[*i].ward_dist(&(clusters[*j]),md.d);
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
    fprintf(stderr, "P delta delta %f\n",delta);
    for (int u : current_clusters)
    {
        for (int v : current_clusters)
        {
            if (u == v)
                continue;
            float distance = clusters[u].ward_dist(&(clusters[v]),delta);
            if (distance <= delta)
            {
                P_delta.push_back(Dist(u,v,distance));
                if (distance < md.d)
                {
                    md.d = distance;
                    md.U = u;
                    md.V = v;
                }
            }
        }
    }
    fprintf(stderr, "P_delta size = %d md = (%d %d %f)\n",P_delta.size(), md.U, md.V, md.d);
    return md;
}
void Clusterizer::ClusterDendrogramm::make(int n)
{
    std::list<Dist> P_delta;
    float delta;
    Dist min = get_P_delta(n,current_clusters,P_delta,delta);
    for (int i=1; i<size; i++)
    {
        if (P_delta.empty())
            min = get_P_delta(n,current_clusters,P_delta,delta);
        else if (min.U < 0)
        {
            for (Dist &d : P_delta)
            {
                if (d.U == d.V)
                    fprintf(stderr,"error in P_delta %d\n",d.U);
                if (d.d < min.d)
                {
                    min.U = d.U;
                    min.V = d.V;
                    min.d = d.d;
                }
            }
            //fprintf(stderr, "min dist (%d %d %f)\n",min.U, min.V, min.d);
        }
        if (min.d > 0.999)
        {
            break;
            //makes no sense to merge clusters with maximum distance between them.
        }
        current_clusters.remove(min.U);
        current_clusters.remove(min.V);
        clusters.push_back(Cluster(&(clusters[min.U]),&(clusters[min.V])));
        int W = clusters.size() - 1;
        auto dit = P_delta.begin();
        while (dit != P_delta.end())
        {
            Dist &d = *dit;
            //fprintf(stderr, "%d %d\n",dit->U, dit->V);
            if (d.U == min.U || d.V == min.U || d.U == min.V || d.V == min.V)
                dit = P_delta.erase(dit);
            else
                dit++; 
        }
        for (int S : current_clusters)
        {
            float d = clusters[W].ward_dist(&(clusters[S]),delta);
            if (d<delta)
            {
                //fprintf(stderr,"new D(%d %d %f)\n",W, S, d);
                P_delta.push_back(Dist(W,S,d));
            }
        }
        current_clusters.push_back(W);
        //fprintf(stderr,"%d %d --> %d dist = %f\n", min.U, min.V, W, min.d);
        min = Dist(-1,-1,1000);
        int sum = 0;
        for (int S : current_clusters)
        {
            //fprintf(stderr,"cluster %d size = %d\n",S,clusters[S].size);
            sum += clusters[S].size;
        }
        //fprintf(stderr,"sum = %d %d \n",sum, size);
    }
    int sum = 0;
    for (int S : current_clusters)
    {
        fprintf(stderr,"cluster %d size = %d\n",S,clusters[S].size);
        sum += clusters[S].size;
    }
    fprintf(stderr,"sum = %d %d \n",sum, size);
}