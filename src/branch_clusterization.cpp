#include "branch_clusterization.h"
#include "billboard_cloud.h"
#include "generated_tree.h"
#include <set>
using namespace glm;
#define DEBUG 0
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
float partial_dist(std::vector<int> &jc, std::vector<float> &matches,  std::vector<float> &weights)
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
        fprintf(stderr,"size = %d joints count by level [", matches.size());
    for (int i=0;i<joint_counts.size();i++)
    {
        joint_counts[i] += bwd2.joint_counts[i];
        fprintf(stderr,"%d ",joint_counts[i]);
    }
        fprintf(stderr,"]\n");
    part_answer.exact = match_joints(b1,b2,matches,joint_counts,min,max);
    part_answer.from = partial_dist(joint_counts,matches,weights);
    part_answer.to = part_answer.exact ? part_answer.from : 1;
    #if DEBUG
    fprintf(stderr," level = %d size = %d dist %f/1 %d/ (%d + %d) matches\n",b1->level,b1->base_seg_n,cur_dist, distances.size(),b1->joints.size() ,
     b2->joints.size());
    #endif
    
    fprintf(stderr," full dist [%f %f] exact = %d\n",part_answer.from, part_answer.to, part_answer.exact);
    
    return part_answer;
}
bool Clusterizer::set_branches(Tree &t, int layer)
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
        ddt.create(i);
        if (branchHeap.branches.size() > 2)
        {
            int i = 0;
            int gistos[20][11]; 
            for (int j1 = 0; j1<20;j1++)
            {
                for (int j2=0;j2<11;j2++)
                {
                    gistos[j1][j2] = 0;
                }
            }
            for (BranchWithData &b : branches)
            {
                for (BranchWithData &b2 : branches)
                {
                    if (b.b == b2.b)
                        continue;
                        fprintf(stderr,"r(%d, %d) = ",0,i);
                        for (int k=1;k<20;k++)
                        {
                            delta = 0.05*k;
                            Answer a = dist(b,b2);
                            gistos[k][(int)(10*a.from)] ++;
                            fprintf(stderr,"%f ", a.from);
                        }
                        i++;
                        fprintf(stderr,"\n");
                }
            }
            fprintf(stderr,"gistogramms\n");
            for (int j1 = 0; j1<20;j1++)
            {
                for (int j2=0;j2<11;j2++)
                {
                    fprintf(stderr,"%d ", gistos[j1][j2]);
                }
                fprintf(stderr,"\n");
            }
        }
        else
        {
           #if DEBUG
              fprintf(stderr,"too few branches for clusterization\n");
           #endif
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
Clusterizer::Answer Clusterizer::cluster_dist_min(Cluster &c1, Cluster &c2, float min = 1.0, float max = 0.0)
{
    float cur_min = 1.0;
    for (BranchWithData &bwd1 : c1.branches)
    {
        for (BranchWithData &bwd2 : c2.branches)
        {
            Answer a = ddt.get(bwd1.pos,bwd2.pos);
            if (a.from > cur_min)
                continue;
            else if (a.exact)
                cur_min = a.from;
            else
            {
                a = dist(bwd1, bwd2, cur_min, 0.0f);
                ddt.set(bwd1.pos,bwd2.pos,a);
                if (a.from<cur_min)
                    cur_min = a.from;
            }
            
        }
    }
}