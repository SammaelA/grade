
#include "generated_tree.h"
#include "../tinyEngine/utility.h"
#include "../clustering/clustering.h"
#include "../texture_manager.h"
#include "../visualizer.h"
#include "../distribution.h"
#include <math.h>
#include <algorithm>
#include "../body.h"
#include <chrono>
#include "../tinyEngine/save_utils/saver.h"
#include "../impostor.h"
#include "../terrain.h"
#include "../field_2d.h"
#include "../grove_generation_utils.h"
#include "../synthetic_trees_generator.h"

namespace mygen
{
#define PI 3.14159265f
int seg_count = 0;
int j_count = 0;
int b_count = 0;
float sum_feed[10];
float count_feed[10];
int ids_counter = 0;
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
        Leaf *l = curTree->leaves->new_leaf();
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
    Branch *nb = curTree->branchHeaps[b->level + 1]->new_branch();
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
    float base_angle = curParams().base_angle();
    glm::vec3 dir = curParams().dir_conserv() * prev_dir + curParams().spread() * N + curParams().phototrop() * M + curParams().gravitrop() * up;
    
    glm::vec3 perp = glm::normalize(dir - prev_dir*glm::dot(prev_dir,dir));
    glm::vec3 base_dir = cos(base_angle)*prev_dir + sin(base_angle)*perp;
    dir += curParams().base_angle_q()*base_dir;

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
            s.mults[i] = MAX(s.mults[i], 3);
            s.mults[i] = 1;
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
    float dist_power = curParams().dist_power();
    float dist_mul = curParams().dist_mul();
    float dist_base = curParams().base_dist();
    float dist_fine = MAX(dist_mul*powf(base->joints.back().distance - dist_base,dist_power), 0);
    float feed = base->joints.back().light - dist_fine;

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
    float remove_chance = std::exp(-average_feed/curParams().branch_removal());
    if (b->size && (b->base_seg_n == 0) && (b->level > 0) && dice(remove_chance,1))
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
    float l = j.childBranches.empty() ? voxels->get_occlusion(j.pos) : 1e6;
    float base_light = curParams().base_light();
    float light_to_feed_pow = curParams().base_light_pow();
    l = base_light*powf(base_light/(base_light + l),light_to_feed_pow);
    if (l < 1)
        l = 1;
    j.light = l;
    for (auto b : j.childBranches)
    {
        b->distance = j.distance;
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
    float dist = b->distance;
    glm::vec3 pr_pos = b->joints.front().pos;
    for (auto &j : b->joints)
    {
        dist += glm::length(pr_pos - j.pos);
        pr_pos = j.pos;
        j.distance = dist;
        l += calc_light(j);
    }
    b->light = l;
    return l;
}
#define MAX_GEN_JOINTS 128
float weights[128];
void TreeGenerator::recalculate_thickness(Branch *b)
{
    bool pr = (b->level == 0) && (curTree->iter % 10 == 0);
    int i = b->joints.size() - 1;
    weights[i + 1] = 0.5;
    auto rev_it = b->joints.rbegin();
    //rev_it++;
    while (rev_it != b->joints.rend())
    {
        weights[i] = weights[i + 1] + MAX(rev_it->light,0);
        for (auto br : rev_it->childBranches)
        {
            weights[i] += MAX(br->light,0);
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
        float r = sum_r * powf((weights[i+1] + j_it->light) / weights[i - 1], 1 / curParams().r_split_save_pow());
        float min_r = 0.05*curParams().base_r();
        s_it->rel_r_begin = MAX(r,min_r);
        s_prev->rel_r_end = s_it->rel_r_begin;
        for (auto br : j_it->childBranches)
        {
            r = s_it->rel_r_begin * powf(br->light / weights[i - 1], 1 / curParams().r_split_save_pow());
            curParams.set_state(br->level);
            min_r = 0.05*curParams().base_r();
            br->base_r = MAX(r,min_r);
            br->base_r *= curParams().base_r();
            curParams.set_state(b->level);
        }
        if (s_prev->rel_r_begin > 3)
        {
            logerr("%f i weight[i-1] [i] sum_r j_it->light rel_r_begin %d %f %f %f %f %f",
            (weights[i] + j_it->light) / weights[i - 1], i,
                    weights[i-1],weights[i],sum_r,j_it->light,s_it->rel_r_begin);
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
uint last_id = 0;
void TreeGenerator::plant_tree(Tree &t, ParameterSetWrapper params)
{
    for (int i = 0; i < 10; i++)
    {
        sum_feed[i] = 0;
        count_feed[i] = 0;
    }
    t.id = last_id;
    last_id++;
    t.voxels = voxels;
    t.params = params;
    if (t.leaves)
        delete t.leaves;
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
    root->type_id = t.type_id;

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

    curTree = &t;
    curParams = params;
    t.root = root;
}
void TreeGenerator::grow_tree(Tree &t)
{
    curParams = t.params;
    curTree = &t;
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
        root->distance = 0;
        feed = calc_light(root);
        root->light += 1000 + 100*sqrt(feed);
        calc_size(root);
        root->size = 0.01;
        grow_branch(root, feed);
        t.iter++;

        if (!(t.iter % 10))
        {
            debugl(2,"tree is growing  %d/%d iteration %d leaves\n", t.iter, (int)curParams().growth_iterations(), 
                  (int)(t.leaves->leaves.size()));
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
bool TreeGenerator::is_branch_productive(Branch *b)
{
    return false && (b->level < curParams().max_depth() - 1) && (b->max_seg_count > b->segments.size());
}

void TreeGenerator::create_grove(Tree *trees, int count)
{
    float r = sqrt(count);
    glm::vec3 vox_center = glm::vec3(0, 100, 0) + curGgd.pos;
    glm::vec3 vox_size = curGgd.size;
    ParameterSetWrapper params = trees[0].params;
    params.set_state(params().max_depth() - 1);
    float single_voxel_size = params().seg_len_mult()/ params().light_precision();
    
    voxels = new LightVoxelsCube(vox_center, vox_size, params().seg_len_mult(), params().light_precision());

    voxels->add_heightmap(*heightmap);
    for (int i=0;i<curGgd.obstacles.size();i++)
    {
        voxels->add_body(curGgd.obstacles[i]);
        seeder->add_body(curGgd.obstacles[i]);
    }
    voxels->calculte_precise_occlusion_from_bodies();
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
void pack_branch_recursively(::Branch *b, GrovePacked &grove, std::vector<unsigned> &ids, BranchStructure &b_struct, int lvl_from, int lvl_to)
{
    if (b->level > lvl_to)
        return;
    if (b->level >= lvl_from)
    {
        PackedBranch pb;
        b->pack(pb);
        ids.push_back(grove.instancedCatalogue.add(pb, b->level));
        b_struct.childBranches.push_back(BranchStructure(ids.back()));
    }
    for (::Joint &j : b->joints)
    {
        for (::Branch *br : j.childBranches)
            pack_branch_recursively(br, grove, ids, b_struct.childBranches.back(), lvl_from, lvl_to);
    }
}
void TreeGenerator::reset()
{
    if (seeder)
        delete seeder;
    if (voxels)
        delete voxels;
    root = nullptr;
    test = nullptr;
    heightmap = nullptr;
    seg_count = 0;
    j_count = 0;
    b_count = 0;
}
void TreeGenerator::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    reset();
    Tree trees[MAX_TREES];
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    heightmap = &h;
    curGgd = ggd;
    seeder = new Seeder(ggd,10,&h);
    int synts = ggd.synts_count;
    int count = ggd.trees_count;
    TreeStructureParameters default_params;
    for (int i = 0; i < count; i++)
    {
        int k = i % ggd.types.size();
        auto &type = ggd.types[k];
        trees[i] = Tree();
        trees[i].type_id = type.type_id;
        TreeStructureParameters *ptr = dynamic_cast<TreeStructureParameters*>(type.params);
        if (!ptr)
        {
            logerr("mygen::TreeGenerator : type of parameters set and generator mismatch. Default parameters will be used.");
            ptr = &default_params;
        }
        trees[i].params = ParameterSetWrapper(*ptr,ptr->max_depth() + 1);
        trees[i].wood = type.wood;
        trees[i].leaf = type.leaf;
    }
    for (int i = 0; i < count + synts; i++)
    {
        int k = i % ggd.types.size();
        auto &type = ggd.types[k];
        trees_external[i] = ::Tree();
        trees_external[i].type = &(ggd.types[k]);
    }

    create_grove(trees, count);

    for (int i = 0; i < count; i++)
    {
        post_process(ggd, trees[i]);
    }

    convert(trees,trees_external,count);

    delete voxels;
    voxels = nullptr;
    delete seeder;
    seeder = nullptr;
    debugl(10,"created %d joints %d branches totally\n",j_count, b_count);
}
void down_stripe(std::vector<float> &res, float start, float end, int count, float pw,float sigma)
{
    res.push_back(start);
    if (count == 1)
        return;
    Normal *gen = distibutionGenerator.get_normal(1,sigma);
    float p, dp = 1/(float)(count - 1);
    for (int i=1;i<count;i++)
    {
        float rel = gen->get()*pow(1-i*dp,pw);
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
        //logerr("no deforms");
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
Tree::Tree():
wood(textureManager.empty()),
leaf(textureManager.empty())
{
    leaves = nullptr;
}
Tree::~Tree()
{
    if (leaves)
    {
        delete leaves;
        leaves = nullptr;
    }
    for (int i=0;i<models.size();i++)
        delete models[i];
    for (int i=0;i<billboardClouds.size();i++)
        delete billboardClouds[i];
    for (int i=0;i<branchHeaps.size();i++)
    {
        delete branchHeaps[i];
    }
}
    void TreeGenerator::convert(mygen::Tree *src, ::Tree *dst, int count)
    {
        for (int i=0;i<count;i++)
        {
            convert(src[i],dst[i]);
        }
    }
    void TreeGenerator::convert(mygen::Tree &src, ::Tree &dst)
    {
        dst.clear();
        dst.id = src.id;
        dst.leaves = new ::LeafHeap();
        dst.pos = src.pos + glm::vec3(0,-000,0);
        for (int i = 0;i<src.branchHeaps.size();i++)
        {
            dst.branchHeaps.push_back(new ::BranchHeap());
        }
        if (src.root)
        {
            dst.root = dst.branchHeaps[src.root->level]->new_branch();
            convert(src, dst,*src.root,*dst.root);
            glm::mat4 trans_matrix = glm::translate(glm::mat4(1.0f),glm::vec3(0,-000,0));
            dst.root->transform(trans_matrix);
            dst.valid = true;
        }
    }
    void TreeGenerator::convert(mygen::Tree &src_tree, ::Tree &dst_tree, mygen::Branch &src, ::Branch &dst)
    {
        dst.center_par = src.center_par;
        dst.center_self = src.center_self;
        dst.id = src_tree.id;
        dst.self_id = ids_counter;
        ids_counter++;
        dst.level = src.level;
        dst.plane_coef = src.plane_coef;
        dst.type_id = src.type_id;
        for (Joint &j : src.joints)
        {
            dst.joints.emplace_back();
            convert(src_tree,dst_tree,j,dst.joints.back());
        }
        for (Segment &s : src.segments)
        {
            dst.segments.emplace_back();
            convert(src_tree,dst_tree,s,dst.segments.back());
        }
    }
    void TreeGenerator::convert(mygen::Tree &src_tree, ::Tree &dst_tree, mygen::Joint &src, ::Joint &dst)
    {
        dst.pos = src.pos;
        if (src.leaf && !src.leaf->dead)
        {
            ::Leaf *l = dst_tree.leaves->new_leaf();
            convert(src_tree,dst_tree,*src.leaf,*l);
            dst.leaf = l;
        }
        for (auto *b : src.childBranches)
        {
            if (b->dead)
                continue;
            ::Branch *br = dst_tree.branchHeaps[b->level]->new_branch();
            convert(src_tree,dst_tree,*b,*br);
            dst.childBranches.push_back(br);
        }
    }
    void TreeGenerator::convert(mygen::Tree &src_tree, ::Tree &dst_tree, mygen::Segment &src, ::Segment &dst)
    {
        dst.begin = src.begin;
        dst.end = src.end;
        dst.mults = src.mults;
        dst.rel_r_begin =src.rel_r_begin;
        dst.rel_r_end = src.rel_r_end;
    }
    void TreeGenerator::convert(mygen::Tree &src_tree, ::Tree &dst_tree, mygen::Leaf &src, ::Leaf &dst)
    {
        dst.pos = src.pos;
        dst.type = 0;
        dst.edges = src.edges;
    }
}
//namespace mygen