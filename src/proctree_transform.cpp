#include "proctree.h"
#include "tree.h"
#include "terrain.h"

float scl = 1;
void br_transform(Proctree::Tree &src_tree, ::Tree &dst_tree, Proctree::Branch *src, ::Branch *dst, int level)
{
    Proctree::fvec3 from_p = src->mParent ? src->mParent->mHead : Proctree::fvec3{0,0,0};
    Proctree::fvec3 to_p = src->mHead;
    
    glm::vec3 from = glm::vec3(from_p.x,from_p.y,from_p.z);
    glm::vec3 to = glm::vec3(to_p.x,to_p.y,to_p.z);
    glm::vec3 dir = glm::normalize(to - from);

    dst->center_par = from;
    dst->center_self = to;
    dst->level = level;
    dst->type_id = dst_tree.type->type_id;
    glm::vec3 N = glm::vec3(1,0,0);
    if (abs(glm::dot(N,dir)) > 1 - 1e-3)
        N = glm::vec3(0,1,0);
    glm::vec3 plane_n = glm::normalize(glm::cross(N,dir));
    float d = -glm::dot(plane_n, to);
    dst->plane_coef = glm::vec4(plane_n.x, plane_n.y, plane_n.z, d);
    if (!src->mChild0 && !src->mChild1)
    {
        dst->segments.emplace_back();
        dst->segments.back().begin = from;
        dst->segments.back().end = to;
        dst->segments.back().rel_r_begin = scl*src->mRadius;
        dst->segments.back().rel_r_end = 0.1*scl*src->mRadius;

        dst->joints.emplace_back();
        dst->joints.back().pos = from;

        dst->joints.emplace_back();
        dst->joints.back().pos = to;
    }

    else if (!src->mChild0 || !src->mChild1)
    {
        Proctree::Branch *chb_src = src->mChild0 ? src->mChild0 : src->mChild1;

        dst->segments.emplace_back();
        dst->segments.back().begin = from;
        dst->segments.back().end = to;
        dst->segments.back().rel_r_begin = scl*src->mRadius;
        dst->segments.back().rel_r_end = scl*src->mRadius;

        dst->joints.emplace_back();
        dst->joints.back().pos = from;

        dst->joints.emplace_back();
        dst->joints.back().pos = to;
        ::Branch *chb = dst_tree.branchHeaps[level]->new_branch();
        br_transform(src_tree,dst_tree,chb_src,chb,level);
        dst->joints.back().childBranches.push_back(chb);
    }
    else
    {
        auto to1 = from + 0.95f*(to - from);
        auto to2 = from + 1.05f*(to - from);

        dst->segments.emplace_back();
        dst->segments.back().begin = from;
        dst->segments.back().end = to1;
        dst->segments.back().rel_r_begin = scl*src->mRadius;
        dst->segments.back().rel_r_end = scl*src->mRadius;

        dst->segments.emplace_back();
        dst->segments.back().begin = to1;
        dst->segments.back().end = to2;
        dst->segments.back().rel_r_begin = scl*src->mRadius;
        dst->segments.back().rel_r_end = scl*src->mRadius;
        dst->joints.emplace_back();
        dst->joints.back().pos = from;

        dst->joints.emplace_back();
        dst->joints.back().pos = to1;
        ::Branch *chb1 = dst_tree.branchHeaps[level + 1]->new_branch();
        src->mHead = Proctree::fvec3{to1.x,to1.y,to1.z};
        br_transform(src_tree,dst_tree,src->mChild0,chb1,level+1);
        dst->joints.back().childBranches.push_back(chb1);

        dst->joints.emplace_back();
        dst->joints.back().pos = to2;
        ::Branch *chb2 = dst_tree.branchHeaps[level + 1]->new_branch();
        src->mHead = Proctree::fvec3{to2.x,to2.y,to2.z};
        br_transform(src_tree,dst_tree,src->mChild1,chb2,level+1);
        dst->joints.back().childBranches.push_back(chb2);
    }

}
void Proctree::transform(Proctree::Tree &src, ::Tree &dst, glm::vec3 pos, glm::vec3 scale)
{
    scl = scale.x;
    int segments = src.mProperties.mSegments;
    
    dst.leaves = new LeafHeap();
    for (int i=0;i<src.mProperties.mLevels + 25;i++)
    {
        dst.branchHeaps.push_back(new BranchHeap());
    }
    dst.pos = glm::vec3(0,0,0);
    dst.root = dst.branchHeaps[0]->new_branch();
    br_transform(src,dst,src.mRoot,dst.root,0);
    
    glm::mat4 trans = glm::scale(glm::translate(glm::mat4(1.0f),pos),scale);
    dst.root->transform(trans);
    
    delete src.mRoot;
	src.mRoot = 0;
}
void Proctree::create_grove(GroveGenerationData ggd, ::Tree *trees, Heightmap &h)
{

    int synts = ggd.synts_count;
    int count = ggd.trees_count;
    for (int i = 0; i < count + synts; i++)
    {
        int k = i % ggd.types.size();
        auto &type = ggd.types[k];
        trees[i] = ::Tree();
        trees[i].type = &(ggd.types[k]);
    }   
    for (int i = 0; i < count; i++)
    {
        Proctree::Tree tree;
        basic_use(tree);
        glm::vec3 pos = glm::vec3(100*i,0,0);
        pos.y = h.get_height(pos);
        transform(tree,trees[i],pos,glm::vec3(40,40,40));
    }
}