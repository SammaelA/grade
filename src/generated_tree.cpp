
#include "generated_tree.h"
#include "tinyEngine/utility.h"
#include "branch_clusterization.h"
#include "texture_manager.h"
#include "visualizer.h"
#include "distribution.h"
#include <math.h>
#include <algorithm>
#include "body.h"
#include <chrono>
#include "tinyEngine/save_utils/saver.h"
#include "impostor.h"
#include "terrain.h"
#include "field_2d.h"
#include "grove_generation_utils.h"
#include "synthetic_trees_generator.h"

#define PI 3.14159265f
int seg_count = 0;
int j_count = 0;
int b_count = 0;
int iter = 0;
float sum_feed[10];
float count_feed[10];
bool dice(float x, float base)
{
    float f = urand();
    return f < (x / base);
}
bool dice_min(float x, float min)
{
    if (x < 1e-8)
        x = 1e-8;
    return dice(min, x);
}
glm::vec3 TreeGenerator::rand_dir()
{
    glm::vec3 dir;
    dir.x = urand(-1,1);
    dir.y = urand(-1,1);
    dir.z = urand(-1,1);
    return glm::normalize(dir);
}
glm::vec3 rand_planar_dir()
{
    glm::vec3 dir;
    dir.x = urand(-1,1);
    dir.y = urand(-1,1);
    dir.z = urand(-1,1);
    return glm::normalize(dir);
}
void TreeGenerator::new_joint(Branch *b, Joint &j)
{
    if (b->level == curParams().max_depth() - 1)
    {
        j.type = j.LEAF;
        Leaf *l = curTree.leaves->new_leaf();
        l->pos = j.pos;
        glm::vec3 rd1 = rand_dir();
        glm::vec3 rd2 = rand_dir();

        glm::vec3 a = j.pos + curParams().leaf_size_mult() * rd1 + 0.5f * curParams().leaf_size_mult() * rd2;
        glm::vec3 b = j.pos + 0.5f * curParams().leaf_size_mult() * rd2;
        glm::vec3 c = j.pos - 0.5f * curParams().leaf_size_mult() * rd2;
        glm::vec3 d = j.pos + curParams().leaf_size_mult() * rd1 - 0.5f * curParams().leaf_size_mult() * rd2;
        l->edges.push_back(a);
        l->edges.push_back(b);
        l->edges.push_back(c);
        l->edges.push_back(d);

        j.leaf = l;
    }
    else if (b->joints.empty())
    {
        j.type = j.MIDDLE;
    }
    else
    {
        float b_ch = powf((b->base_seg_n + (float)b->joints.size() / (b->base_seg_n + b->max_seg_count)),
                          curParams().branching_power());
        b_ch = curParams().min_branching_chance() + b_ch * (curParams().max_branching_chance() - curParams().min_branching_chance());
        if (dice(b_ch, 1.0))
        {
            j.type = j.FORK;
            j.max_branching = 1 + (curParams().max_branching() - 1) * floor(urand());
        }
        else
            j.type = j.MIDDLE;
    }
    voxels->set_occluder(j.pos, powf(curParams().max_depth() - b->level,2));
    calc_light(j);
    b->joints.push_back(j);
    j_count++;
}
void TreeGenerator::new_branch(Branch *b, Joint &j, Segment &s, glm::vec3 &M, bool from_end)
{
    Branch *nb = curTree.branchHeaps[b->level + 1]->new_branch();
    nb->base_seg_n = from_end ? b->segments.size() : 0;
    nb->level = b->level + 1;
    nb->max_seg_count = curParams().max_segments();
    nb->base_r = curParams().base_r();
    nb->type_id = b->type_id;
    Joint nj;
    nj.pos = j.pos;
    new_joint(nb, nj);
    glm::vec3 prev_dir = glm::normalize(s.end - s.begin);
    glm::vec3 rnd = rand_planar_dir();
    glm::vec3 N = glm::normalize(rnd - dot(rnd, prev_dir) * rnd); //Normal Vector
    glm::vec3 up = glm::vec3(0, 1, 0);

    glm::vec3 dir = curParams().dir_conserv() * prev_dir + curParams().spread() * N + curParams().phototrop() * M + curParams().gravitrop() * up;
    dir = from_end ? prev_dir : glm::normalize(dir);

    new_segment2(nb, dir, j.pos);
    calc_light(nb);
    calc_size(nb);
    nb->light += j.light;
    j.light = 0;

    glm::vec3 plane_n = glm::normalize(glm::cross(N,dir));
    float d = -dot(plane_n, j.pos);
    nb->plane_coef = glm::vec4(plane_n.x, plane_n.y, plane_n.z, d);

    j.childBranches.push_back(nb);
    test = nb;
    b_count++;
}
int branch_try_schedule(int k)
{
    return k <= 3 ? 1 : powf(2,(k-3)/2) + 1;
}
void TreeGenerator::try_new_branch(Branch *b, Joint &j, Segment &s, bool from_end)
{
    curParams.set_state(b->level + 1);
    float feed = j.light;
    int bs = j.childBranches.size();
    if ((b->level < curParams().max_depth() - 1) && 
        ((j.type == j.FORK && &j != &(b->joints.back()) && 
        (bs < j.max_branching)) || from_end))
    {
        j.iters_to_next_branch_try--;
        if (j.iters_to_next_branch_try <= 0)
        {
            j.tries_from_last_grown_branch++;
            if (dice(feed, exp(bs) * curParams().base_branch_feed()))
            {
                float len = curParams().max_segments()*curParams().seg_len_mult();
                float occ = 0.0;
                glm::vec3 M;
                if (b->level == 0)
                    M = voxels->get_dir_to_bright_place_cone(j.pos,0.33*len,256, &occ);
                else if (b->level == 1)
                    M = voxels->get_dir_to_bright_place_cone(j.pos,0.33*len, 64, &occ);
                else if (b->level == 2)
                    M = voxels->get_dir_to_bright_place_ext(j.pos, 3, &occ);
                else 
                    M = voxels->get_dir_to_bright_place_ext(j.pos, 1, &occ);
                if (dice(1, occ * curParams().branch_grow_decrease_q()) && (occ < 100000))
                {
                    new_branch(b, j, s, M, from_end);
                    j.tries_from_last_grown_branch = 0;
                }
            }
            j.iters_to_next_branch_try = branch_try_schedule(j.tries_from_last_grown_branch);
        }
    }
    curParams.set_state(b->level);
}
inline float lfun(float v1, float v2, float a, float b)
{
    return (v1 / v2) * (b - a) + a;
}
void TreeGenerator::set_seg_r(Branch *base, Segment &s)
{
    float b_st = curParams().base_r();
    if (base->level + 1 < curParams().max_depth())
        curParams.set_state(base->level + 1);
    float b_end = curParams().base_r();
    curParams.set_state(base->level);
    int n = base->segments.size();
    s.rel_r_end = lfun(curParams().max_segments() - n - 1, curParams().max_segments(), b_end, b_st);
    s.rel_r_begin = lfun(curParams().max_segments() - n - 1, curParams().max_segments(), b_end, b_st);
}
void TreeGenerator::new_segment2(Branch *base, glm::vec3 &dir, glm::vec3 &pos)
{
    Segment s;
    float len = curParams().seg_len_mult();
    s.begin = pos;
    s.end = s.begin + len * dir;
    set_seg_r(base, s);
    if (base->level < curParams().r_deformation_levels())
    {
        int deforms = curParams().r_deformation_points();
        if (base->segments.empty() || base->segments.back().mults.size() != deforms)
            s.mults = std::vector<float>(deforms,1);
        else
        {
            s.mults = base->segments.back().mults;
        }
        
        for (int i=0;i<deforms;i++)
        {
            float k = curParams().r_deformation_power();
            k = k > 0 ? 1 + k : 1/(1-k);
            s.mults[i] *= k;
        }
    }
    Joint j;
    j.pos = s.end;
    base->segments.push_back(s);
    new_joint(base, j);
}
void TreeGenerator::new_segment(Branch *base, glm::vec3 &M)
{
    Joint &last = base->joints.back();
    Segment &last_s = base->segments.back();

    float sp = curParams().seg_spread();
    float bn = base->level * curParams().seg_bend() * powf((float)base->segments.size() / base->max_seg_count, curParams().seg_bend_pow());
    glm::vec3 prev_dir = glm::normalize(last_s.end - last_s.begin);
    glm::vec3 rnd = rand_dir();
    glm::vec3 N = glm::normalize(glm::cross(prev_dir, rnd)); //Normal Vector
    glm::vec3 dir = glm::normalize(glm::mix(prev_dir, N, sp));
    glm::vec3 up = glm::vec3(0, 1, 0);

    dir = curParams().seg_dir_conserv() * prev_dir + curParams().seg_spread() * N + curParams().seg_phototrop() * M + curParams().seg_gravitrop() * up - bn * up;
    dir = glm::normalize(dir);
    float q = 0;
    dir = get_optimal_segment_growth_direction(q,base);
    if (glm::length(glm::cross(dir,glm::vec3(base->plane_coef))) < 0.001)
        dir = glm::normalize(dir + glm::vec3(0.05,0,0));
    if (glm::length(glm::cross(dir,glm::vec3(base->plane_coef))) < 0.001)
        dir = glm::normalize(dir + glm::vec3(0.0,0.05,0));
    new_segment2(base, dir, last.pos);
}
int segment_try_schedule(int k)
{
    return k <= 3 ? 1 : 1.5*k - 4;
}
void TreeGenerator::try_new_segment(Branch *base)
{
    float feed = base->joints.back().light;
    if (base->segments.size() < base->max_seg_count)
    {
        base->iters_to_next_segment_try--;
        if (base->iters_to_next_segment_try <= 0)
        {
            base->tries_from_last_grown_segment++;
            if (dice(feed, curParams().base_seg_feed()))
            {
                float occ = 0.0;
                sum_feed[base->level] += feed;
                count_feed[base->level] += 1;
                glm::vec3 M = voxels->get_dir_to_bright_place_ext(base->segments.back().end, 1, &occ);
                if (dice(1, occ * curParams().segment_grow_decrease_q()) && (occ < 100000))
                {
                    new_segment(base, M);
                    base->tries_from_last_grown_segment = 0;
                }
            }
            base->iters_to_next_segment_try = segment_try_schedule(base->tries_from_last_grown_segment);
        }
    }
}
bool comp(float *a, float *b)
{
    return *a > *b;
}
void TreeGenerator::distribute_feed(Branch *b)
{
    std::vector<float *> feeds;
    float *top = &(b->joints.back().light);
    //feeds.push_back(top);
    float min_light = 1.05;
    for (auto &j : b->joints)
    {
        if (j.childBranches.empty() && j.light > min_light)
            feeds.push_back(&(j.light));
        for (auto ch_b : j.childBranches)
        {
            if (ch_b->light > min_light)
                feeds.push_back(&(ch_b->light));
        }
    }

    std::sort(feeds.begin(), feeds.end(), comp);

    float w = 1.0;
    float w_sum = 0.0;
    float w_decay = 1 - curParams().feed_distribution_d_weight();
    for (auto f : feeds)
    {
        (*f) = 1 / w;
        w_sum += 1 / w;
        w++;
        if (1 / w < curParams().feed_distribution_min_weight())
            w = 1 / curParams().feed_distribution_min_weight();
    }
    float top_bonus = is_branch_productive(b) ? w_sum * curParams().top_growth_bonus() / (1 - curParams().top_growth_bonus()) : 0;
    w_sum += top_bonus;
    *top += top_bonus;
    float mult = b->light / w_sum;
    for (auto f : feeds)
    {
        (*f) *= mult;
    }
}
void TreeGenerator::remove_branch(Branch *b)
{
    b->dead = true;
    bool save_trunk = false;
    Joint j1, j2;
    Segment s;
    if (b->joints.size() >= 2 && b->segments.size() >= 1)
    {
        save_trunk = true;
        auto jj = b->joints.begin();
        j1 = *jj;
        j2 = *(jj++);
        s = b->segments.front();
    }
    for (auto j : b->joints)
    {
        if (j.leaf)
            j.leaf->dead = true;
        for (auto br : j.childBranches)
            remove_branch(br);
    }
    b->joints.clear();
    b->segments.clear();
    if (save_trunk)
    {
        b->joints.push_back(j1);
        b->joints.push_back(j2);
        b->segments.push_back(s);
    }
}
void TreeGenerator::grow_branch(Branch *b, float feed)
{
    if (b->dead)
        return;
    curParams.set_state(b->level);
    float average_feed = b->light / (b->size + 0.001);
    if (b->size && (b->base_seg_n == 0) && !dice(average_feed, curParams().branch_removal()))
    {
        remove_branch(b);
        return;
    }
    distribute_feed(b);

    auto j = b->joints.begin();
    j++;
    auto s = b->segments.begin();
    while ((j != b->joints.end()) && (s != b->segments.end()))
    {
        bool from_end = (b->joints.size() == b->max_seg_count + 1) && (std::next(j) == b->joints.end()) && (j->childBranches.size() == 0);
        try_new_branch(b, *j, *s, from_end);
        j++;
        s++;
    }
    try_new_segment(b);
    recalculate_thickness(b);
    for (auto &j : b->joints)
    {
        for (auto br : j.childBranches)
        {
            grow_branch(br, 1);
        }
    }
}
float TreeGenerator::calc_light(Joint &j)
{
    float l = j.childBranches.empty() ? voxels->get_occlusion(j.pos) : 0;
    l = 50 - l;
    if (l < 1)
        l = 1;
    j.light = l;
    for (auto b : j.childBranches)
    {
        l += calc_light(b);
    }
    return l;
}
float TreeGenerator::calc_size(Joint &j)
{
    float sz = 0;
    for (auto b : j.childBranches)
    {
        sz += calc_size(b);
    }
    return sz;
}
float TreeGenerator::calc_size(Branch *b)
{
    float sz = b->segments.size() * powf(2, curParams().max_depth() - b->level);
    for (auto &j : b->joints)
    {
        sz += calc_size(j);
    }
    b->size = sz;
    return sz;
}
float TreeGenerator::calc_light(Branch *b)
{
    if (b->dead)
        return 0;
    float l = 0.0;
    for (auto &j : b->joints)
    {
        l += calc_light(j);
    }
    b->light = l;
    return l;
}
void TreeGenerator::recalculate_thickness(Branch *b)
{
    float *weights = new float[b->joints.size() + 1];
    int i = b->joints.size() - 1;
    weights[i + 1] = 0.5;
    auto rev_it = b->joints.rbegin();
    while (rev_it != b->joints.rend())
    {
        weights[i] = weights[i + 1] + rev_it->light;
        for (auto br : rev_it->childBranches)
        {
            weights[i] += br->light;
        }
        i--;
        rev_it++;
    }
    i = 1;
    auto s_it = b->segments.begin();
    auto s_prev = b->segments.begin();
    auto j_it = b->joints.begin();
    j_it++;
    s_it++;
    s_prev->rel_r_begin = b->base_r;
    while (s_it != b->segments.end())
    {
        float sum_r = s_prev->rel_r_begin;
        s_it->rel_r_begin = sum_r * powf((weights[i] + j_it->light) / weights[i - 1], 1 / curParams().r_split_save_pow());
        s_prev->rel_r_end = s_it->rel_r_begin;
        for (auto br : j_it->childBranches)
        {
            br->base_r = s_it->rel_r_begin * powf(br->light / weights[i - 1], 1 / curParams().r_split_save_pow());
            curParams.set_state(br->level);
            br->base_r *= curParams().base_r();
            curParams.set_state(b->level);
        }
        s_prev++;
        s_it++;
        i++;
        j_it++;
    }
    s_prev->rel_r_end = s_prev->rel_r_begin * powf(weights[i] / weights[i - 1], 1 / curParams().r_split_save_pow());
    if (!b->joints.back().childBranches.empty())
    {
        Branch *cont = b->joints.back().childBranches.front();
        curParams.set_state(cont->level);
        s_prev->rel_r_end = MIN(s_prev->rel_r_begin, s_prev->rel_r_end*curParams().base_r());
        cont->base_r = s_prev->rel_r_end;
         curParams.set_state(b->level);
    }
    delete[](weights);
}
glm::vec3 TreeGenerator::get_optimal_branch_growth_direction(float &quality, Branch *base, Joint &j, Segment &s, bool from_end)
{

}
glm::vec3 TreeGenerator::get_optimal_segment_growth_direction(float &quality, Branch *base)
{
    LightVoxelsCube *field = nullptr;
    glm::vec3 prev_dir = glm::normalize(base->segments.back().end - base->segments.back().begin);
    float bn = base->level * curParams().seg_bend() * powf((float)base->segments.size() / base->max_seg_count, curParams().seg_bend_pow());
    
    calc_quality_field(field,base->joints.back().pos,glm::vec3(5*(curParams().max_depth() - base->level)),prev_dir,base->plane_coef,
                      curParams().seg_dir_conserv(),
                      curParams().seg_plane_conserv(),
                      curParams().seg_dir_random(),
                      curParams().seg_spread(),
                      curParams().seg_gravitrop() - bn,
                      curParams().seg_phototrop());
    
    glm::ivec3 vox_sizes = field->get_vox_sizes();
    std::vector<glm::vec3> best_pos;
    float best = 1e9;
    for (int i = -vox_sizes.x; i <= vox_sizes.x; i++)
    {
        for (int j = -vox_sizes.y; j <= vox_sizes.y; j++)
        {
            for (int k = -vox_sizes.z; k <= vox_sizes.z; k++)
            {
                float q = field->get_occlusion_voxel(glm::ivec3(i,j,k));

                if (abs(q - best) < 0.01)
                {
                    best_pos.push_back(glm::vec3(i,j,k));
                }
                else if (q < best)
                {
                    best = q;
                    best_pos.clear();
                    best_pos.push_back(glm::vec3(i,j,k));
                }
            }
        }
    }
    int p = urandi(0,best_pos.size());
    quality = best;
    if (field)
        delete field;
    return glm::normalize(best_pos[p]);
}
void TreeGenerator::calc_quality_field(LightVoxelsCube *&field, glm::vec3 pos, glm::vec3 sizes, glm::vec3 prev_dir, glm::vec4 plane,
                                       float dir_cons, float plane_cons, float rnd, float spread, float up, float to_light)
{
    field = new LightVoxelsCube(voxels,pos,sizes);
    glm::ivec3 vox_sizes = field->get_vox_sizes();
    float voxel_size = field->get_voxel_size();
    uint64_t w = 0,x = 0, s = (uint64_t)(1e4*pos.x + 1e6*pos.y + 1e8*pos.z);
    for (int i = -vox_sizes.x; i <= vox_sizes.x; i++)
    {
        for (int j = -vox_sizes.y; j <= vox_sizes.y; j++)
        {
            for (int k = -vox_sizes.z; k <= vox_sizes.z; k++)
            {
                glm::ivec3 vx = glm::ivec3(i,j,k);
                glm::vec3 dir = glm::normalize(glm::vec3(i,j,k));
                float len = voxel_size*sqrt(i*i + j*j + k*k);

                float rnd_q = srand(s,x,w);
                float dir_cons_q = 1 - dot(dir,prev_dir);
                float plane_cons_q = abs(dot(dir,glm::vec3(plane.x,plane.y,plane.z)));
                float up_q =  1 - dot(dir, glm::vec3(0,1,0));
                float spread_q = dot(dir,prev_dir);
                float to_light_q = field->get_occlusion_voxel(vx);

                float quality = rnd*rnd_q + dir_cons*dir_cons_q + plane_cons*plane_cons_q + up*up_q + 
                                spread*spread_q + to_light*to_light_q;
                field->replace_occluder_voxel(vx,quality);
            }
        }
    }
}
LightVoxelsCube *TreeGenerator::create_light_voxels_cube(ParameterSetWrapper params, glm::vec3 pos)
{
    float max_size = 0.0;
    float wm = 0.7;
    float hm = 0.55;
    for (int i = 0; i < params().max_depth(); i++)
    {
        params.set_state(i);
        max_size += params().max_segments() * params().seg_len_mult();
    }
    params.set_state(params().max_depth() - 1);
    return new LightVoxelsCube(pos + glm::vec3(0, 0.5 * max_size, 0), glm::vec3(wm * max_size, hm * max_size, wm * max_size),
                               params().seg_len_mult(), params().light_precision());
}
void TreeGenerator::plant_tree(Tree &t, ParameterSetWrapper params)
{
    for (int i = 0; i < 10; i++)
    {
        sum_feed[i] = 0;
        count_feed[i] = 0;
    }
    t.voxels = voxels;
    glm::vec3 sun_dir(-1, -1, -1);
    t.params = params;
    LeafHeap *lh = new LeafHeap();
    t.leaves = lh;
    for (int i = 0; i < params().max_depth(); i++)
    {
        BranchHeap *bh = new BranchHeap();
        t.branchHeaps.push_back(bh);
    }
    curParams.set_state(0);
    root = t.branchHeaps[0]->new_branch();
    root->level = 0;
    root->max_seg_count = curParams().max_segments();
    root->base_seg_n = 0;
    root->base_r = curParams().base_r();
    root->type_id = t.type ? t.type->type_id : 0;

    Joint j1, j2;
    Segment ts;
    j1.pos = t.pos;
    new_joint(root, j1);
    ts.begin = j1.pos;
    ts.end = j1.pos + glm::vec3(0, 3*curParams().base_r(), 0);
    j2.pos = ts.end;
    ts.rel_r_begin = 1;
    ts.rel_r_end = 1;

    new_joint(root, j2);
    root->segments.push_back(ts);
    root->plane_coef = glm::normalize(glm::vec4(urand(0,1),0,urand(0,1),0));
    root->plane_coef.w = -glm::dot(glm::vec3(root->plane_coef),root->joints.front().pos);

    curTree = t;
    curParams = params;
    t.root = root;
}
void TreeGenerator::grow_tree(Tree &t)
{
    curParams = t.params;
    curTree = t;
    root = t.root;
    if (t.branchHeaps.size() > 1 && !t.branchHeaps[1]->branches.empty())
        test = &(t.branchHeaps[1]->branches.front());
    test = nullptr;
    voxels = t.voxels;
    for (int i = 0; i < 10; i++)
    {
        sum_feed[i] = 0;
        count_feed[i] = 0;
    }
    if (t.iter < curParams().growth_iterations() && root)
    {
        float feed = 1;
        seg_count = 0;
        feed = calc_light(root);
        root->light += 1000 + 100*sqrt(feed);
        calc_size(root);
        root->size = 0.01;
        grow_branch(root, feed);
        t.iter++;

        if (!(t.iter % 10))
        {
            debugl(2,"tree is growing  %d/%d iteration %d leaves\n", t.iter, curParams().growth_iterations(), t.leaves->leaves.size());
            debugl(2, "sum feed %f\n", feed);
            debugl(2, "average feed distribution:");
            for (int j = 0; j < 10; j++)
            {
                if (count_feed[j] < 0.1)
                    debugl(2, " nan");
                else
                {
                    debugl(2, " %f", sum_feed[j] / count_feed[j]);
                }
            }
            debugl(2, "\n");
            voxels->print_average_occlusion();
        }
    }
}
bool TreeGenerator::tree_to_model(Tree &t, bool leaves, DebugVisualizer &debug)
{
    int sz = 1024;
    for (int level = 0; level < 3; level++)
    {
        Model *m = new Model();
        for (int i = 0; i < level; i++)
        {
            for (auto &branch : t.branchHeaps[i]->branches)
            {
                debug.branch_to_model(branch, m, leaves);
            }
        }
        if (leaves)
        {
            for (auto &leaf : t.leaves->leaves)
            {
                debug.leaf_to_model(leaf, m);
            }
        }
        m->scale = glm::vec3(t.params().scale());
        BillboardCloudRaw *cloud = new BillboardCloudRaw(sz, sz, curGgd.types);
        cloud->prepare(t, level, level);
        t.billboardClouds.push_back(cloud);
        sz = sz * 2;
        t.models.push_back(m);
    }
    BillboardCloudRaw *cloud = new BillboardCloudRaw(1, 1, curGgd.types);
    Model *m = new Model();
    for (int i = 0; i < t.branchHeaps.size(); i++)
    {
        for (auto &branch : t.branchHeaps[i]->branches)
        {
            debug.branch_to_model(branch, m, leaves);
        }
    }

    Model *m2 = new Model();
    if (leaves)
    {
        for (auto &leaf : t.leaves->leaves)
        {
            debug.leaf_to_model(leaf, m2);
        }
    }
    t.billboardClouds.push_back(cloud);
    t.models.push_back(m);
    t.models.push_back(m2);
    return true;
}
void LeafHeap::clear_removed()
{
}
void BranchHeap::clear_removed()
{
}
bool TreeGenerator::is_branch_productive(Branch *b)
{
    return false && (b->level < curParams().max_depth() - 1) && (b->max_seg_count > b->segments.size());
}
void TreeGenerator::create_tree(Tree &t, ParameterSetWrapper params, DebugVisualizer &debug)
{
    voxels = create_light_voxels_cube(params, t.pos);
    plant_tree(t, params);
    while (t.iter < params().growth_iterations())
    {
        grow_tree(t);
    }
    //Clusterizer cl;
    //cl.set_branches(t,2, debug);
    debug.set_params(t.params());
    //debug.add_branch(t.root,glm::vec3(1,1,1),glm::vec3(0,100,0),3);
    tree_to_model(t, false, debug);
}
void TreeGenerator::create_grove(Tree *trees, int count, DebugVisualizer &debug)
{
    float r = sqrt(count);
    glm::vec3 vox_center = glm::vec3(0, 100, 0) + curGgd.pos;
    glm::vec3 vox_size = curGgd.size;
    ParameterSetWrapper params = trees[0].params;
    params.set_state(params().max_depth() - 1);
    Box b = Box(glm::vec3(12,0,-50),glm::vec3(100,0,0),glm::vec3(0,100,0),glm::vec3(0,0,100));
    Cylinder el = Cylinder(glm::vec3(250,50,-20),glm::vec3(70,0,0),glm::vec3(0,50,0),glm::vec3(0,0,50));
    float single_voxel_size = params().seg_len_mult()/ params().light_precision();
    
    voxels = new LightVoxelsCube(vox_center, vox_size, params().seg_len_mult(), params().light_precision());

    voxels->add_heightmap(*heightmap);
    for (int i=0;i<curGgd.obstacles.size();i++)
    {
        voxels->add_body(curGgd.obstacles[i]);
        seeder->add_body(curGgd.obstacles[i]);
    }
    const int growth_step = 10;
    int trees_planted = 0;  
    for (int j = 0; j < params().growth_iterations(); j++)
    {
        if (j % growth_step == 0)
        {
            seeder->recalcuate_shadows(trees,trees_planted);

            int cnt = MAX((count - trees_planted < 3) ? count - trees_planted : (count - trees_planted)/2,0);
            std::vector<Seed> seeds;
            seeder->choose_places_for_seeds(cnt,seeds);
            for (Seed &seed : seeds)
            {
                glm::vec3 pos = glm::vec3(seed.pos.x,0,seed.pos.y);
                pos.y = heightmap->get_height(pos);
                
                trees[trees_planted].pos = pos;
                plant_tree(trees[trees_planted], trees[trees_planted].params);
                trees_planted++;
            }
        }
        
        for (int i = 0; i < trees_planted; i++)
        {
            if (trees[i].iter < params().growth_iterations())
            {
                grow_tree(trees[i]);
            }
        }
    }
    for (int j = 0; j < params().growth_iterations(); j++)
    {
        for (int i = 0; i < trees_planted; i++)
        {
            if (trees[i].iter < params().growth_iterations())
            {
                grow_tree(trees[i]);
            }
        }
    }
}
void pack_branch_recursively(Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, BranchStructure &b_struct, int lvl_from, int lvl_to)
{
    if (b->dead)
        return;
    if (b->level > lvl_to)
        return;
    if (b->level >= lvl_from)
    {
        PackedBranch pb;
        b->pack(pb);
        ids.push_back(grove.instancedCatalogue.add(pb, b->level));
        b_struct.childBranches.push_back(BranchStructure(ids.back()));
    }
    for (Joint &j : b->joints)
    {
        for (Branch *br : j.childBranches)
            pack_branch_recursively(br, grove, ids, b_struct.childBranches.back(), lvl_from, lvl_to);
    }
}
void pack_cluster(ClusterData &cluster, GrovePacked &grove, std::vector<BranchStructure> &instanced_structures, int lvl_from, int lvl_to)
{
    instanced_structures.push_back(BranchStructure());
    grove.instancedBranches.push_back(InstancedBranch());
    grove.instancedBranches.back().IDA = cluster.IDA;
    std::vector<unsigned> &ids = grove.instancedBranches.back().branches;
    Branch *base = cluster.base;
    pack_branch_recursively(base, grove, ids,instanced_structures.back(), lvl_from, lvl_to);
    if (instanced_structures.back().childBranches.size() == 1)
        instanced_structures.back() = instanced_structures.back().childBranches[0];

    for (int i=0;i<cluster.ACDA.originals.size();i++)//leave marks on branch to construct tree structure in future
    {
        cluster.ACDA.originals[i]->base_seg_n = instanced_structures.size() - 1;//cluster id
        cluster.ACDA.originals[i]->max_seg_count = - i - 100;//i is number in cluster
    }
    grove.instancedBranches.back().bbox = BillboardCloudRaw::get_minimal_bbox(base);
}
void pack_tree(Tree &t, GrovePacked &grove, int up_to_level)
{
    for (int i = 0; i <= up_to_level; i++)
    {
        for (Branch &branch : t.branchHeaps[i]->branches)
        {
            PackedBranch b;
            branch.pack(b);
            branch.max_seg_count = grove.uniqueCatalogue.add(b, branch.level);
        }
    }
}
void pack_structure(Branch *rt, GrovePacked &grove, BranchStructure &str, std::vector<BranchStructure> &instanced_structures)
{
    str.pos = rt->max_seg_count;
    for (Joint &j : rt->joints)
    {
        for (Branch *br : j.childBranches)
        {
            if (br->max_seg_count < 0)//it is a mark made by pack_cluster
            {
                unsigned transform_n = - (br->max_seg_count + 100);
                unsigned instance_n = br->base_seg_n;
                BranchStructure bs = instanced_structures[instance_n];
                glm::mat4 tr = grove.instancedBranches[instance_n].IDA.transforms[transform_n];
                str.childBranchesInstanced.push_back(std::pair<glm::mat4,BranchStructure>(tr,bs));
            }
            else
            {
                BranchStructure bs;
                pack_structure(br,grove,bs,instanced_structures);
                str.childBranches.push_back(bs);
            }
            
        }
    }
}

void transform_according_to_root(ClusterData &cluster)
{
    InstanceDataArrays IDA = cluster.IDA;
    Branch *base = cluster.base;
    if (cluster.ACDA.originals.size() != IDA.transforms.size())
    {
        logerr("Transform according to root failed: corrupted clusters data.");
        return;
    }
    for (int i = 0; i< cluster.ACDA.originals.size(); i++)
    {
        if (base->joints.empty() || cluster.ACDA.originals[i]->joints.empty())
            continue;
        std::vector<glm::vec3> replaced_joints;
        for (Joint &j : base->joints)
            replaced_joints.push_back(glm::vec3(IDA.transforms[i]*glm::vec4(j.pos,1)));
        for (Joint &j : cluster.ACDA.originals[i]->joints)
        {
            float d_min = 1e9;
            glm::vec3 n_min = glm::vec3(0);
            for (glm::vec3 &np : replaced_joints)
            {
                float d = glm::length(np - j.pos);
                if (d<d_min)
                {
                    d_min = d;
                    n_min = np;
                }
            }
            glm::vec3 sh = n_min - j.pos;
            //logerr("shifting after trunk clusterization %f %f %f",sh.x,sh.y,sh.z);
            for (Branch *ch_b : j.childBranches)
            {
                glm::mat4 tr = glm::translate(glm::mat4(1.0f),sh);
                ch_b->transform(tr);
            }
        }
    }
}
void TreeGenerator::create_grove(GroveGenerationData ggd, GrovePacked &grove, DebugVisualizer &debug, Tree *trees,
                                 Heightmap *h, bool visualize_voxels)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    heightmap = h;
    curGgd = ggd;
    seeder = new Seeder(ggd,10,h);
    int synts = ggd.synts_count;
    int count = ggd.trees_count;
    for (int i = 0; i < count + synts; i++)
    {
        int k = i % ggd.types.size();
        auto &type = ggd.types[k];
        trees[i] = Tree();
        trees[i].type = &(ggd.types[k]);
        trees[i].params = ParameterSetWrapper(type.params,type.params.max_depth() + 1);
        trees[i].wood = type.wood;
        trees[i].leaf = type.leaf;
    }
    grove.center = glm::vec3(0,0,0);
    grove.ggd_name = ggd.name;

    create_grove(trees, count, debug);

    for (int i = 0; i < count; i++)
    {
        post_process(ggd, trees[i]);
    }

    Clusterizer tr_cl;
    ClusterizationParams tr_cp;
    tr_cp.weights = std::vector<float>{1,0,0,0.0,0.0};
    tr_cp.ignore_structure_level = 1;
    tr_cp.delta = 0.05;
    tr_cp.light_importance = 0;
    tr_cp.different_types_tolerance = true;
    tr_cp.r_weights = std::vector<float>{0.5,0,0,0.0,0.0}; 
    tr_cp.max_individual_dist = 0.1;
    tr_cp.bwd_rotations = 4;
    tr_cl.set_clusterization_params(tr_cp);
    tr_cl.set_branches(trees, count, 0, voxels);

    std::vector<ClusterData> trunks_clusters, branches_clusters, full_tree_clusters;
    tr_cl.clusterize(trunks_clusters);
    for (ClusterData &cd : trunks_clusters)
    {
        transform_according_to_root(cd);
    }
    Clusterizer cl;
    ClusterizationParams cp;
    cp.weights = std::vector<float>{5000,800,40,0.0,0.0};
    cp.ignore_structure_level = 1;
    cp.delta = 0.3;
    cp.max_individual_dist = 0.9;
    cp.bwd_rotations = 4;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    cl.set_clusterization_params(cp);
    cl.set_branches(trees, count, 1, voxels);

    cl.clusterize(branches_clusters);
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    for (int i = 0;i<branches_clusters.size();i++)
    {
        for (Branch *br : branches_clusters[i].ACDA.originals)
        {
            br->base_seg_n = i;
        }
    }

    SyntheticTreeGenerator stg = SyntheticTreeGenerator(*seeder, trunks_clusters, branches_clusters, curGgd);
    stg.generate(trees + count, synts, voxels);
    
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    
    Clusterizer cl2;
    cp.ignore_structure_level = 1;
    cp.light_importance = 0.8;
    cp.different_types_tolerance = false;
    cl2.set_clusterization_params(cp);
    cl2.set_branches(trees, count + synts, 0, voxels);
    cl2.clusterize(full_tree_clusters);

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    grove.clouds.push_back(BillboardCloudData());
    grove.clouds.push_back(BillboardCloudData());//empty 'zero' data

    ImpostorBaker *ib = new ImpostorBaker(2048, 2048, curGgd.types);
    grove.impostors.push_back(ImpostorsData());
    ib->prepare_all_grove(trees[0], ggd, 0, full_tree_clusters, &grove.impostors.back());
    grove.impostors.push_back(ImpostorsData());
    ImpostorBaker *ib2 = new ImpostorBaker(BillboardCloudRaw::Quality::MEDIUM,0,full_tree_clusters,curGgd.types,&grove.impostors.back());

    
    grove.clouds.push_back(BillboardCloudData());//main cloud 1
    BillboardCloudRaw *cloud1 = new BillboardCloudRaw(BillboardCloudRaw::Quality::HIGH, 1,
                                                      branches_clusters,curGgd.types,&grove.clouds.back());
    
    grove.clouds.push_back(BillboardCloudData());//main cloud 2
    BillboardCloudRaw *cloud2 = new BillboardCloudRaw(BillboardCloudRaw::Quality::LOW, 2,
                                                      branches_clusters,curGgd.types,&grove.clouds.back());
    std::vector<BranchStructure> instanced_structures;

    for (ClusterData &cd : trunks_clusters)
    {
        pack_cluster(cd, grove, instanced_structures, 0, 0);
    }
    for (ClusterData &cd : branches_clusters)
    {
        pack_cluster(cd, grove, instanced_structures, 1, 1000);
    }
    for (int i = 0; i < count; i++)
    {
        grove.roots.push_back(BranchStructure());
        pack_structure(trees[i].root,grove,grove.roots.back(),instanced_structures);
    }
    if (visualize_voxels)
    {
        debug.visualize_light_voxels(voxels,ggd.pos-ggd.size + glm::vec3(0,ggd.size.y,0),2.0f*ggd.size - glm::vec3(0,ggd.size.y,0),glm::vec3(2),0.8,0.1);
    }
    delete(voxels);
    for (int i = 0; i < count; i++)
    {
        trees[i].voxels = nullptr;
    }

    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
    std::cerr << "Generation took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "[ms]" << std::endl;
    std::cerr << "Main clusterization took " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "[ms]" << std::endl;
    std::cerr << "Syntetic trees generation took " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << "[ms]" << std::endl;
    std::cerr << "Secondary clusterization took " << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count() << "[ms]" << std::endl;
    std::cerr << "Finishing took " << std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count() << "[ms]" << std::endl;
    logerr("created %d joints %d branches totally",j_count, b_count);
}
void down_stripe(std::vector<float> &res, float start, float end, int count, float pw,float sigma)
{
    res.push_back(start);
    if (count == 1)
        return;
    Normal gen = Normal(1,sigma);
    float p, dp = 1/(float)(count - 1);
    for (int i=1;i<count;i++)
    {
        float rel = gen.get()*pow(1-i*dp,pw);
        rel = rel*(start - end) + end;
        if (rel < end)
            rel = end;
        if (rel > res.back())
            rel = res.back();
        res.push_back(rel);
    }
}
void TreeGenerator::deform_root(Branch *b)
{
    if (!b)
        return;
    if (b->segments.empty())
        return;
    int s = b->segments.back().mults.size();
    if (s == 0)
    {
        logerr("no deforms");
        return;
    }
    int stripes = 3;
    for (int i=0;i<stripes;i++)
    {
        int p = s*i/stripes;
        std::vector<float> res;
        down_stripe(res,2.5,1,b->segments.size(),1.75,0.05);
        int k = 0;
        for (auto &seg : b->segments)
        {
            if (urandi(0,b->segments.size()/2) == 0)
                p = (p + 1) % s;
            if (k >= res.size())
                break;
            if (seg.mults.size() < s)
            {
                seg.mults = std::vector<float>(s,1);
            }
            seg.mults[p] *= res[k];
            float dr = 1 + 0.4*(res[k] - 1);
            float dr2 = 1 + 0.4*(dr - 1);
            seg.mults[(p - 1 + s) % s] *= dr;
            seg.mults[(p + 1 + s) % s] *= dr;
            seg.mults[(p - 2 + s) % s] *= dr2;
            seg.mults[(p + 2 + s) % s] *= dr2;
            k++;
        }
    }
}
void TreeGenerator::post_process(GroveGenerationData ggd, Tree &t)
{
    deform_root(t.root);
    set_branches_centers(ggd, t, 1);
}
glm::vec3 b_center_self(Branch *b)
{
    glm::vec3 c;
    for (Joint &j : b->joints)
        c += j.pos;
    c = c /(float)b->joints.size();
    return c;
}
void b_center_rec(Branch *b, glm::vec3 center_par, int up_to_level)
{
    b->center_par = center_par;
    if (b->level <= up_to_level)
        b->center_self = 0.5f*(b->joints.front().pos + b->joints.back().pos);
    else
        b->center_self = b->center_par;
    for (Joint &j : b->joints)
    {
        for (Branch *br : j.childBranches)
        {
            b_center_rec(br,b->center_self,up_to_level);
        }
    }
}
void TreeGenerator::set_branches_centers(GroveGenerationData ggd, Tree &t, int up_to_level)
{
    if (!t.root)
        return;
    b_center_rec(t.root,ggd.pos,up_to_level);
}