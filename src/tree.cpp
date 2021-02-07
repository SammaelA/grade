#include "tree.h"
#include "texture_manager.h"
#include "billboard_cloud.h"
Tree::Tree():
wood(textureManager.empty()),
leaf(textureManager.empty())
{
    leaves = nullptr;
}
Tree::~Tree()
{
    if (leaves)
        leaves->leaves.clear();
    for (int i=0;i<models.size();i++)
        delete models[i];
    //for (int i=0;i<branchHeaps.size();i++)
    //    delete branchHeaps[i];
    for (int i=0;i<billboardClouds.size();i++)
        delete billboardClouds[i];
}
void Branch::deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap)
{
    level = b->level;
    id = b->id;
    type_id = b->type_id;
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
        for (Branch *br : j.childBranches)
        {
            Branch *nb = heap.new_branch();
            nb->deep_copy(br, heap, leaf_heap);
            nj.childBranches.push_back(nb);
        }
        joints.push_back(nj);
    }
}
void Branch::transform(glm::mat4 &trans_matrix)
{
    for (Segment &s : segments)
    {
        s.begin = trans_matrix * glm::vec4(s.begin, 1.0f);
        s.end = trans_matrix * glm::vec4(s.end, 1.0f);
    }
    for (Joint &j : joints)
    {
        j.pos = trans_matrix * glm::vec4(j.pos, 1.0f);
        if (j.leaf && !(j.leaf->dead))
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
    if (dead || joints.size() <= 1 || (joints.size() != segments.size() + 1))
        return;

    auto jit = joints.begin();
    auto sit = segments.begin();

    branch.joints.clear();
    branch.leaves.clear();
    branch.level = level;
    branch.type_id = type_id;
    branch.joints.push_back(PackedJoint(sit->begin, sit->rel_r_begin));
    jit++;
    bool has_deforms = false;
    while (jit != joints.end())
    {
        PackedLeaf l;
        if (jit->leaf && !(jit->leaf->dead) && (jit->leaf->edges.size() >= 3))
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
TreeTypeData::TreeTypeData(int id, TreeStructureParameters _params, std::string wood_tex_name, std::string leaf_tex_name):
wood(textureManager.get(wood_tex_name)),
leaf(textureManager.get(leaf_tex_name))
{
    type_id = id;
    params = _params;
}
float Branch::get_r_mult(float phi, std::vector<float> &mults)
{
    if (mults.empty())
        return 1;
    else 
    {
        int b = (int)(phi/(2*PI)*mults.size());
        float q = phi/(2*PI)*mults.size() - b;

        return (1-q)*mults[b % mults.size()] + q*mults[(b + 1) % mults.size()];
    }
}