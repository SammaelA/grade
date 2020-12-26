#include "branch_clusterization.h"
#include "billboard_cloud.h"
#include "generated_tree.h"
#include <set>
using namespace glm;
#define DEBUG 0
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
    if (b.length()<0.01)
        b = vec3(0,1,0);
    else
        b = normalize(b);
    c = cross(a,b);

    bbox = BillboardCloud::get_bbox(branch,a,b,c);
    return true;
}
float partial_dist(std::vector<int> &jc, std::vector<float> &matches,  const std::vector<float> &weights)
{
    float num=0.0;
    float denom = 0.0;
    for (int i=0;i<matches.size();i++)
    {
        num += weights[i]*(2*matches[i]);
        denom +=weights[i]*jc[i];
    }
    return denom > 0.001 ? 1 - num/denom : 0;
}
bool Clusterizer::match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc, float min, float max)
{
    int sz1 = j1->childBranches.size();
    int sz2 = j2->childBranches.size();
    if (sz1 == 0 || sz2 == 0)
        return true;
    else if (sz1 == 1 && sz2 == 1)
    {
        return match_joints(j1->childBranches.front(), j2->childBranches.front(), matches, jc, min, max);
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
bool Clusterizer::match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, float min, float max)
{
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
            it++;
        }
    }
    //after it only correct matches remained here
    int matches_count = distances.size();
    matches[b1->level] += matches_count;
    if (partial_dist(jc,matches,weights) > min)
        return false;
    if (b1->level == matches.size() - 1)//this is the last branch level, no need to iterate over child branches
        return true;
    it = distances.begin();
    while (it != distances.end())
    {
        if (!match_child_branches(it->j1,it->j2,matches,jc,min,max))
        {
            //if child branches matching made early exit, it means that min distance limit
            //is already passed. No need to go further.
            return false;
        }
        it++;
    }
    return true;

}
Clusterizer::Answer Clusterizer::dist(BranchWithData &bwd1, BranchWithData &bwd2, float min, float max)
{
    Branch *b1 = bwd1.b;
    Branch *b2 = bwd2.b;
    Answer part_answer;
    std::vector<int> joint_counts(bwd1.joint_counts);
    std::vector<float> matches(joint_counts.size());
    //fprintf(stderr,"size = %d joints count by level [", matches.size());
    for (int i=0;i<joint_counts.size();i++)
    {
        joint_counts[i] += bwd2.joint_counts[i];
        //fprintf(stderr,"%d ",joint_counts[i]);
    }
    //fprintf(stderr,"]\n");
    part_answer.exact = match_joints(b1,b2,matches,joint_counts,min,max);
    part_answer.from = partial_dist(joint_counts,matches,weights);
    part_answer.to = part_answer.exact ? part_answer.from : 1;
    #if DEBUG
    fprintf(stderr," level = %d size = %d dist %f/1 %d/ (%d + %d) matches\n",b1->level,b1->base_seg_n,cur_dist, distances.size(),b1->joints.size() ,
     b2->joints.size());
    #endif
    
    //fprintf(stderr," full dist [%f %f] exact = %d\n",part_answer.from, part_answer.to, part_answer.exact);
    
    return part_answer;
}
bool Clusterizer::set_branches(Tree &t, int layer, DebugVisualizer &debug)
{
    if (layer<0 || layer>t.branchHeaps.size() || t.branchHeaps[layer]->branches.size() == 0)
    {
        #if DEBUG
            fprintf(stderr,"no branches for clusterization\n");
        #endif
        return false;
    }
    else
    {
        int i = 0;
        for (Branch &b : t.branchHeaps[layer]->branches)
        {
            //if (b.dead)
            //    continue;
            BillboardCloud::BBox bbox;
            Branch *nb = branchHeap.new_branch();
            nb->deep_copy(&b,branchHeap);
            //nb = &b;
            if (dedicated_bbox(nb,bbox)) 
            {
            
                mat4 rot_inv(vec4(bbox.a,0),vec4(bbox.b,0),vec4(bbox.c,0),vec4(0,0,0,1));
                mat4 rot = inverse(rot_inv);
                vec3 base_joint_pos = rot*vec4(b.joints.front().pos,1.0f);
                mat4 transl = translate(mat4(1.0f),-1.0f*base_joint_pos);
                mat4 SC = scale(mat4(1.0f), 0.01f*bbox.sizes);
                mat4 SC_inv = inverse(SC);
                rot = SC_inv*transl*rot;
                transform_branch(nb,rot);
                branches.push_back(BranchWithData(nb,t.params.max_depth(),i));
                i++;
                
            }
        }
        ClusterDendrogramm Ddg;
        Ddg.make_base_clusters(branches);
        Ddg.make(20);
        std::vector<Branch *>branches;
        int k = 0;
        for (int S : Ddg.current_clusters)
        {
            Ddg.clusters[S].to_branch_data(branches);
            for (int i=0;i<branches.size();i++)
            {
                debug.add_branch(branches[i],vec3(1,1,1),vec3(k,0,0),layer);
                debug.add_branch(branches[i],vec3(1,1,1),vec3(k,100,0),layer);
                debug.add_branch(branches[i],vec3(1,1,1),vec3(k,100,0),layer+1);
            }
            branches.clear();
            k+=100;
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
Clusterizer::ClusterDendrogramm::Dist 
Clusterizer::ClusterDendrogramm::get_P_delta(int n,std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta)
{
    fprintf(stderr, "get P delta\n");
    Dist md(-1,-1,1000);
    int k = n > current_clusters.size() ? current_clusters.size() : n;
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
            float distance = clusters[*i].ward_dist(&(clusters[*j]));
            if (distance < md.d)
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
            float distance = clusters[u].ward_dist(&(clusters[v]));
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
    for (Dist &d : P_delta)
    {
        //fprintf(stderr, "(%d %d %f)\n",d.U, d.V, d.d);
    }
    return md;
}
void Clusterizer::ClusterDendrogramm::make(int n)
{
    std::list<Dist> P_delta;
    float delta;
    Dist min = get_P_delta(n,current_clusters,P_delta,delta);
    for (int i=1;i<size - n;i++)
    {
        if (P_delta.empty())
            min = get_P_delta(n,current_clusters,P_delta,delta);
        else if (min.U < 0)
        {
            for (Dist &d : P_delta)
            {
                if (d.d < min.d)
                    min.U = d.U;
                    min.V = d.V;
                    min.d = d.d;
            }
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
            float d = clusters[W].ward_dist(&(clusters[S]));
            if (d<0.0001)
                return;
            if (d<delta)
            {
                fprintf(stderr,"new D(%d %d %f)\n",W, S, d);
                P_delta.push_back(Dist(W,S,d));
            }
        }
        current_clusters.push_back(W);
        fprintf(stderr,"%d %d --> %d dist = %f\n", min.U, min.V, W, min.d);
        min = Dist(-1,-1,1000);
    }
    for (int S : current_clusters)
    {
        fprintf(stderr,"cluster %d size = %d\n",S,clusters[S].size);
    }
}