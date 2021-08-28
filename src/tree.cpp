#include "tree.h"
#include "texture_manager.h"
#include "billboard_cloud.h"
#include "spline.h"

void Branch::norecursive_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap)
{
    
    level = b->level;
    id = b->id;
    type_id = b->type_id;
    mark_A = b->mark_A;
    mark_B = b->mark_B;
    dead = b->dead;
    plane_coef = b->plane_coef;
    center_par = b->center_par;
    center_self = b->center_self;

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

        if (leaf_heap && j.leaf)
        {
            Leaf *l = leaf_heap->new_leaf();
            *l = *(j.leaf);
            nj.leaf = l;
        }
        else
        {
            nj.leaf = nullptr; //TODO: it's temporal
        }
        nj.childBranches.clear();

        joints.push_back(nj);
    }
}
void Branch::deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap)
{
    if (!b)
        return;
    level = b->level;
    id = b->id;
    type_id = b->type_id;
    mark_A = b->mark_A;
    mark_B = b->mark_B;
    dead = b->dead;
    plane_coef = b->plane_coef;
    center_par = b->center_par;
    center_self = b->center_self;

    segments.clear();
    joints.clear();
    for (const Segment &s : b->segments)
    {
        Segment seg = s;
        segments.push_back(seg);
    }
    for (const Joint &j : b->joints)
    {
        Joint nj;
        nj.pos = j.pos;
        nj.mark_A = j.mark_A;
        if (leaf_heap && j.leaf)
        {
            Leaf *l = leaf_heap->new_leaf();
            *l = *(j.leaf);
            nj.leaf = l;
        }
        else
        {
            nj.leaf = nullptr; //TODO: it's temporal
        }
        for (Branch *br : j.childBranches)
        {
            Branch *nb = heap.new_branch();
            nb->deep_copy(br, heap, leaf_heap);
            nj.childBranches.push_back(nb);
        }
        joints.push_back(nj);
    }
}
void Branch::transform(glm::mat4 &trans_matrix, float r_transform)
{
    plane_coef.w = 0;
    plane_coef = glm::normalize(trans_matrix * plane_coef);
    plane_coef.w = -dot(glm::vec3(plane_coef),joints.front().pos);
    for (Segment &s : segments)
    {
        s.begin = trans_matrix * glm::vec4(s.begin, 1.0f);
        s.end = trans_matrix * glm::vec4(s.end, 1.0f);
        s.rel_r_begin = r_transform * s.rel_r_begin;
        s.rel_r_end = r_transform * s.rel_r_end;
    }
    for (Joint &j : joints)
    {
        j.pos = trans_matrix * glm::vec4(j.pos, 1.0f);
        if (j.leaf)
        {
            j.leaf->pos = trans_matrix * glm::vec4(j.leaf->pos, 1.0f);
            for (glm::vec3 &vec : j.leaf->edges)
            {
                vec = trans_matrix * glm::vec4(vec, 1.0f);
            }
        }
        for (Branch *br : j.childBranches)
        {
            br->transform(trans_matrix);
        }
    }
}
void Branch::pack(PackedBranch &branch)
{
    if (joints.size() <= 1 || (joints.size() != segments.size() + 1))
        return;

    auto jit = joints.begin();
    auto sit = segments.begin();

    branch.joints.clear();
    branch.leaves.clear();
    branch.level = level;
    branch.type_id = type_id;
    branch.plane_coef = plane_coef;
    branch.joints.push_back(PackedJoint(sit->begin, sit->rel_r_begin));
    jit++;
    bool has_deforms = false;
    while (jit != joints.end())
    {
        PackedLeaf l;
        if (jit->leaf && (jit->leaf->edges.size() >= 3))
        {
            l.edges = jit->leaf->edges;
        }
        branch.joints.push_back(PackedJoint(sit->end, sit->rel_r_end));
        branch.leaves.push_back(l);
        branch.r_mults.push_back(sit->mults);
        has_deforms = has_deforms || !sit->mults.empty();
        jit++;
        sit++;
    }
    if (!has_deforms)
    {
        branch.r_mults.clear();
    }
}
TreeTypeData::TreeTypeData(int id, ParametersSet *_params, std::string wood_tex_name, std::string leaf_tex_name):
wood(textureManager.get(wood_tex_name)),
leaf(textureManager.get(leaf_tex_name))
{
    type_id = id;
    params = _params;
}
TreeTypeData::~TreeTypeData()
{
    //if (params)
    //    delete params;
}
float Branch::get_r_mult(float phi, std::vector<float> &mults)
{
    if (mults.empty())
        return 1;
    else
    {
        std::vector<double> x,y;
        for (int i=0;i<mults.size();i++)
        {
            x.push_back(2*PI*i/mults.size());
            y.push_back(mults[i]);
        }
        x.push_back(2*PI);
        y.push_back(mults[0]);
        auto s = spline(x,y);
        int b = (int)(phi/(2*PI)*mults.size());
        return s[b].get(phi);
    } 
}

void BranchHeap::clear_removed()
{
    auto it = branches.begin();
    logerr("clearing branches: before %d",branches.size());
    while (it != branches.end())
    {
        if (it->dead)
        {
            it = branches.erase(it);
        }
        else
        {
            it++;
        }
    }
    logerr("clearing branches: after %d",branches.size());
}

void LeafHeap::clear_removed()
{
    auto it = leaves.begin();
    logerr("clearing leaves: before %d",leaves.size());
    while (it != leaves.end())
    {
        if (it->edges.empty())
        {
            it = leaves.erase(it);
        }
        else
        {
            it++;
        }
    }
    logerr("clearing leaves: after %d",leaves.size());
}

void Branch::mark_dead()
{
    dead = true;
    for (Joint &j : joints)
    {
        if (j.leaf)
            j.leaf->edges.clear();
        for (Branch *br : j.childBranches)
            br->mark_dead();
    }
}