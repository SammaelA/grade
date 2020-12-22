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
Clusterizer::Answer Clusterizer::dist(Branch *b1, Branch *b2, float min, float max)
{
    if (!b1 || !b2)
    {
        return Answer(true,0,0);
    }
    if (b1->joints.size() == 1)
    {
        if (b2->joints.size() == 1)
            return Answer(true,0,0);
        else
            return Answer(true,1,1);
        
    }
    /*BillboardCloud::BBox bbox1, bbox2;
    if (dedicated_bbox(b1,bbox1) && dedicated_bbox(b2,bbox2)) 
    {

    }
    else
    {
        return Answer(true,0,0);//all broken branches are similar to each other
    }*/
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
    int matches = matched_joints.size();
    cur_dist = 0.5*(b1->joints.size() + b2->joints.size() - matches)/(float)(b1->joints.size() + b2->joints.size());
    Answer part_answer(true,cur_dist,1 - (0.5-cur_dist));
    float subbranches_j_count = b1->base_seg_n + b2->base_seg_n - b1->joints.size() - b2->joints.size();
    subbranches_j_count *= 2.0;
    int matched_j = 0;
    if (subbranches_j_count < 1)
    {
        #if DEBUG
            fprintf(stderr," fast dist %f exact = 1",2*cur_dist);
        #endif
        return Answer(true,2*cur_dist,2*cur_dist);
    }
    it = distances.begin();
    while (it != distances.end())
    {
        if (part_answer.from > min || part_answer.to < max)
        {
            part_answer.exact = false;
            break;
        }
        Answer min_d = Answer(false,1,-1);
        int min_sz = 0;
        int sz1 = it->j1->childBranches.size();
        int sz2 = it->j2->childBranches.size();
        if (sz1 == 0 || sz2 == 0)
        {
            
            if (sz1)
            {
                for (Branch *b1 : it->j1->childBranches)
                {
                    matched_j += b1->base_seg_n;
                    part_answer.from += b1->base_seg_n/subbranches_j_count;
                }
            }   
            else if (sz2)
            {
                for (Branch *b2 : it->j2->childBranches)
                {
                    matched_j += b2->base_seg_n;
                    part_answer.from += b2->base_seg_n/subbranches_j_count;
                }
            }
            it++;
            continue;
        }
        else if (sz1 == 1)
        {
            Branch *f = it->j1->childBranches.front();
            for (Branch *f2 : it->j2->childBranches)
            {
                Answer a = dist(f,f2,(min - cur_dist)*subbranches_j_count/(f->base_seg_n + f2->base_seg_n),
                                (max - cur_dist)*subbranches_j_count/(f->base_seg_n + f2->base_seg_n));
                if (a.from<min_d.from)
                {
                    min_d = a;
                    min_sz = f->base_seg_n + f2->base_seg_n;
                }
            }
            float k = min_sz/subbranches_j_count;
            part_answer.from += k*min_d.from;
            part_answer.to -= k*(1 - min_d.to);
        }  
        else if (sz2 == 1)
        {
            Branch *f = it->j2->childBranches.front();
            for (Branch *f2 : it->j1->childBranches)
            {
                Answer a = dist(f,f2,(min - cur_dist)*subbranches_j_count/(f->base_seg_n + f2->base_seg_n),
                                (max - cur_dist)*subbranches_j_count/(f->base_seg_n + f2->base_seg_n));
                if (a.from<min_d.from)
                {
                    min_d = a;
                    min_sz = f->base_seg_n + f2->base_seg_n;
                }
            }
            float k = min_sz/subbranches_j_count;
            part_answer.from += k*min_d.from;
            part_answer.to -= k*(1 - min_d.to);
        }
        else
        {
            it++;
            continue;
        }
        matched_j += min_sz;
        it++;
    }
    float unmatched_joints = (0.5*subbranches_j_count - matched_j)/subbranches_j_count;
    part_answer.from += unmatched_joints;
    #if DEBUG
    fprintf(stderr," level = %d size = %d dist %f/1 %d/ (%d + %d) matches\n",b1->level,b1->base_seg_n,cur_dist, distances.size(),b1->joints.size() ,
     b2->joints.size());
    #endif
    #if DEBUG
        fprintf(stderr," full dist [%f %f] exact = %d\n",part_answer.from, part_answer.to, part_answer.exact);
    #endif
    
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
                
            }
        }
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
            for (Branch &b : branchHeap.branches)
            {
                if (b.level != layer)
                    continue;
                for (Branch &b2 : branchHeap.branches)
                {
                    if (b.level ==  b2.level)
                    {
                        fprintf(stderr,"r(%d, %d) = ",0,i);
                        for (int k=1;k<20;k++)
                        {
                            delta = 0.05*k;
                            Answer a = dist(&(b),&(b2));
                            gistos[k][(int)(10*a.from)] ++;
                            fprintf(stderr,"%f ", a.from);
                        }
                        i++;
                        fprintf(stderr,"\n");
                    }
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