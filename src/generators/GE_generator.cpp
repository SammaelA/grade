#include "GE_generator.h"
using namespace glm;

int GETreeGenerator::ids = 0;
int GETreeGenerator::t_ids = 0;
int GETreeGenerator::iteration = 0;

void GETreeGenerator::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    for (int i=0;i<ggd.trees_count;i++)
    {
        Tree t;
        t.pos = vec3(50*i,0,0);
        t.pos.y = h.get_height(t.pos);
        GETreeParameters params;
        
        create_tree_internal(t,params);
        convert(t, trees_external[ggd.trees_count - i - 1], ggd);
        t_ids++;
    }
}
void GETreeGenerator::create_tree_internal(Tree &t, GETreeParameters &params)
{
    iteration = 0;
    create_initial_trunk(t, params);
    
    //generate tree
    
    set_levels_rec(t, t.root,params,0);
}
void GETreeGenerator::create_initial_trunk(Tree &t, GETreeParameters &params)
{
    t.root.joints.push_back(Joint(t.pos, 1));
    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,3,0), 0.9));
    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,6,0), 0.7));
    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,9,0), 0.53));
    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,12,0), 0.5));
    t.root.joints[3].childBranches.emplace_back();
    auto &b = t.root.joints[3].childBranches.back();
    vec3 p = t.root.joints[3].pos;
    b.joints.push_back(Joint(p,0.3));
    b.joints.push_back(Joint(p + vec3(1,1,0),0.3));
    b.joints.push_back(Joint(p + vec3(1.3,1.7,0),0.3));
    b.joints.push_back(Joint(p + vec3(1.6,2.5,0),0.3));

    b.joints[2].childBranches.emplace_back();
    auto &b1 = b.joints[2].childBranches.back();
    p = b.joints[2].pos;

    b1.joints.push_back(Joint(p,0.1));
    b1.joints.push_back(Joint(p + vec3(0,0.5,0.5),0.1));
    b1.joints.push_back(Joint(p + vec3(0,1,0.6),0.1));
    b1.joints.push_back(Joint(p + vec3(0.3,1.1,0.7),0.1));
}
void GETreeGenerator::set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level)
{
    b.level = level;
    t.max_depth = level;
    for (Joint &j : b.joints)
    {
        for (Branch &br : j.childBranches)
        {
            if (br.alive)
                set_levels_rec(t, br, params, level + 1);
        }
    }
}
void GETreeGenerator::convert(Tree &src, ::Tree &dst, GroveGenerationData &ggd)
{
    for (int j=0;j<src.max_depth + 1;j++)
    {
        BranchHeap *br = new BranchHeap();
        dst.branchHeaps.push_back(br);
    }

    dst.leaves = new LeafHeap();
    dst.id = t_ids;
    dst.pos = src.pos;
    dst.type = &(ggd.types[0]);
    dst.valid = true;

    dst.root = dst.branchHeaps[0]->new_branch();
    dst.root->type_id = 0;
    dst.root->self_id = ids++;
    dst.root->level = 0;
    dst.root->dead = false;
    dst.root->center_self = src.pos;
    dst.root->center_par = vec3(0,0,0);
    dst.root->plane_coef = vec4(1,0,0,-src.pos.x);
    dst.root->id = dst.id;

    convert(src, dst, src.root, dst.root);
}
void GETreeGenerator::convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst)
{
    for (int i=0;i<b_src.joints.size();i++)
    {
        logerr("conv %d %d",b_src.level, i);
        b_dst->joints.emplace_back();
        b_dst->joints.back().pos = b_src.joints[i].pos;
        if (!b_src.joints[i].leaf.edges.empty())
        {
            b_dst->joints.back().leaf = dst.leaves->new_leaf();
            b_dst->joints.back().leaf->pos = b_src.joints[i].pos;
            b_dst->joints.back().leaf->type = 0;
            b_dst->joints.back().leaf->edges = b_src.joints[i].leaf.edges;
        }
        if (i != 0)
        {
            b_dst->segments.emplace_back();
            b_dst->segments.back().begin = b_src.joints[i-1].pos;
            b_dst->segments.back().rel_r_begin = b_src.joints[i-1].r;
            b_dst->segments.back().end = b_src.joints[i].pos;
            b_dst->segments.back().rel_r_end = b_src.joints[i].r;
        }
        for (Branch &chb : b_src.joints[i].childBranches)
        {
            if (!chb.alive)
                continue;
            ::Branch *nb = dst.branchHeaps[chb.level]->new_branch();
            b_dst->joints.back().childBranches.push_back(nb);
            nb->type_id = 0;
            nb->self_id = ids++;
            nb->level = chb.level;
            nb->dead = false;
            nb->center_self = b_dst->joints.back().pos;
            nb->center_par = b_dst->center_self;
            nb->plane_coef = vec4(1,0,0,-b_dst->joints.back().pos.x);
            nb->id = dst.id;

            convert(src, dst, chb, nb);
        }
    }
}