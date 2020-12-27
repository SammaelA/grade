#include "tree.h"
void Branch::deep_copy(const Branch *b, BranchHeap &heap)
{
    level = b->level;
    base_seg_n = b->base_seg_n;
    max_seg_count = b->max_seg_count;
    light = b->light;
    size = b->size;
    base_r = b->base_r;
    dead = b->dead;
    segments.clear();
    joints.clear();
    for (const Segment &s : b->segments)
    {
        Segment seg = s;
        segments.push_back(seg);
    }
    for (const Joint &j : b->joints)
    {
        Joint nj = j;
        nj.childBranches.clear();
        for (Branch *br : j.childBranches)
        {
            Branch *nb = heap.new_branch();
            nb->deep_copy(br, heap);
            nj.childBranches.push_back(nb);
        }
        joints.push_back(nj);
    }
}
void Branch::transform(glm::mat4 &trans_matrix)
{
    for (Segment &s : segments)
    {
        s.begin = trans_matrix*glm::vec4(s.begin,1.0f);
        s.end = trans_matrix*glm::vec4(s.end,1.0f);
    }
    for (Joint &j : joints)
    {
        j.pos = trans_matrix*glm::vec4(j.pos,1.0f);
        for (Branch *br : j.childBranches)
        {
            br->transform(trans_matrix);
        }
    }
}