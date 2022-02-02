#include "GE_generator.h"
#include <algorithm>
using namespace glm;

std::atomic<int> GETreeGenerator::ids(0);
std::atomic<int> GETreeGenerator::t_ids(0);
int GETreeGenerator::joints_limit = 50000;
GETreeParameters GETreeGenerator::defaultParameters = GETreeParameters();

//TreeTypeData def_ttd = TreeTypeData(-1,&(GETreeGenerator::defaultParameters), "wood","leaf");
bool GETreeGenerator::iterate(LightVoxelsCube &voxels)
{
    bool growing = false;
    for (auto &t : trees)
    {
        GETreeParameters *p_ptr = dynamic_cast<GETreeParameters *>(t.type->params);
        GETreeParameters &params = p_ptr ? *(p_ptr) : defaultParameters;
        if (t.status == TreeStatus::SEED)
        {
            create_initial_trunk(t, params);
            //set_occlusion(t.root, voxels, params, 1);
        }
        if (t.status == TreeStatus::GROWING)
        {
            iteration = t.iteration;

            auto &p = params;
            float A = ((float)p.Xm / p.X0 - 1) * exp(-p.r * iteration);
            float dX_dt = p.Xm * p.r * A / SQR(1 + A);
            int max_growth = round(dX_dt);

            //logerr("iteration %d max growth %d", iteration, max_growth);
            if (iteration >= params.max_iterations || max_growth <= 0)
                t.status == TreeStatus::GROWN;
            else
            {
                std::vector<GrowPoint> growth_points;
                sp_data = SpaceColonizationData();
                sp_data.active = true;
                calc_light(t.root, voxels, params);
                //logerr("light %.1f", t.root.total_light);
                float l0 = t.root.total_light;
                t.root.total_light *= MAX((t.root.total_light + 30)/(t.root.total_light + 0.1),1.5);
                float l1 = CLAMP(t.root.total_light, 0, 15*params.Xm);
                float l2 = CLAMP(t.root.total_light - l1, 0, 30*params.Xm);
                float l3 = CLAMP(t.root.total_light - l2, 0, 60*params.Xm);
                float l4 = CLAMP(t.root.total_light - l2, 0, 1000*params.Xm);
                float l_corr = l1 + 0.4*l2 + 0.16*l3 + 0.0*l4;
                distribute_resource(t.root, params, 1.0);
                prepare_nodes_and_space_colonization(t, t.root, params, growth_points, max_growth);
                sp_data.prepare(voxels);
                //logerr("prepared %d grow points %d sp dots", growth_points.size(), sp_data.positions.size());
                grow_nodes(t, params, growth_points, voxels, max_growth);
                recalculate_radii(t, t.root, params);
                remove_branches(t, t.root, params, voxels);

                t.iteration++;
                growing = true;
            }
        }
    }

    return growing;
}
void GETreeGenerator::plant_tree(glm::vec3 pos, TreeTypeData *type)
{
    GETreeParameters *ps = dynamic_cast<GETreeParameters *>(type->params);
    if (!ps)
    {
        logerr("Tree type %d cannot be generated with GE tree generator", type->type_id);
        type->params = &(GETreeGenerator::defaultParameters);
    }
    trees.emplace_back();
    trees.back().pos = pos;
    trees.back().type = type;
}
void GETreeGenerator::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
    int i = 0;
    for (auto &t : trees)
    {
        vec3 min_pos;
        vec3 max_pos;
        GETreeParameters *p_ptr = dynamic_cast<GETreeParameters *>(t.type->params);
        GETreeParameters &params = p_ptr ? *(p_ptr) : defaultParameters;
        set_levels_rec(t, t.root, params, 0);
        min_pos = vec3(1e9,1e9,1e9);
        max_pos = -min_pos;
        create_leaves(t.root, params, 0, voxels);
        vec3 sz = max_pos - min_pos;
        //logerr("tree size %f %f %f max_sz %f %f %f", sz.x, sz.y, sz.z, params.Xm, 1.5f*params.Xm, params.Xm);
        convert(t, trees_external[trees.size() - i - 1]);
        i++;
    }
    trees.clear();
    iteration = 0;
}
void GETreeGenerator::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    GETreeParameters params;
    ggd.types[0].params = &params;

    ivec3 voxels_sizes = ivec3(2 * ggd.trees_count * params.Xm, 3.5 * params.Xm, 2 * params.Xm);
    LightVoxelsCube voxels = LightVoxelsCube(vec3(0,0,0), voxels_sizes, 0.5 * params.ro);
    for (int i = 0; i < ggd.trees_count; i++)
    {
        vec3 pos = vec3( i % 2 ? 15 * ((i + 1)/2) : -15*(i/2), 0, 0);
        plant_tree(pos, &(ggd.types[0]));
    }
    bool grow = true;
    while (grow)
    {
        grow = iterate(voxels);
    }
    finalize_generation(trees_external, voxels);
}

void GETreeGenerator::create_leaves(Branch &b, GETreeParameters &params, int level_from, LightVoxelsCube &voxels)
{
    for (Joint &j : b.joints)
    {
        //min_pos = min(min_pos, j.pos);
        //max_pos = max(max_pos, j.pos);
        if (j.childBranches.empty() && b.level >= level_from && params.leaves_cnt > 0 &&
            j.r < params.base_r * params.leaves_max_r)
        {
            float f_leaves_cnt = params.leaves_cnt * SQR(1 / (0.5 + voxels.get_occlusion_trilinear(j.pos)));
            int l_cnt = f_leaves_cnt;
            if (self_rand() < (f_leaves_cnt - l_cnt))
                l_cnt++;
            for (int i=0;i<l_cnt;i++)
            {
                glm::vec3 rd1 = normalize(vec3(self_rand(-1, 1), self_rand(-params.leaves_angle_a, params.leaves_angle_a), 
                                            self_rand(-1, 1)));
                glm::vec3 rd2 = normalize(vec3(self_rand(-1, 1), self_rand(-params.leaves_angle_b, params.leaves_angle_b),
                                            self_rand(-1, 1)));
                float sz = params.ro * params.leaf_size_mult;
                glm::vec3 a = j.pos + sz * rd1 + 0.5f * sz * rd2;
                glm::vec3 b = j.pos + 0.5f * sz * rd2;
                glm::vec3 c = j.pos - 0.5f * sz * rd2;
                glm::vec3 d = j.pos + sz * rd1 - 0.5f * sz * rd2;

                j.leaf.pos = j.pos;
                j.leaf.edges.push_back(a);
                j.leaf.edges.push_back(b);
                j.leaf.edges.push_back(c);
                j.leaf.edges.push_back(d);
            }
        }
        for (Branch &br : j.childBranches)
        {
            if (br.alive)
                create_leaves(br, params, level_from, voxels);
        }
    }
}
void GETreeGenerator::create_initial_trunk(Tree &t, GETreeParameters &params)
{
    t.root = Branch();
    t.root.level = 0;

    if (params.root_type == 0)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-1,0), 1, iteration, true));
        for (int i = 0; i < 10; i++)
        {
            float phi = 0.2*PI*(i + self_rand());
            t.root.joints.push_back(Joint(t.pos + glm::vec3(0, 0.1*(i+1)*params.ro, 0), 0.9, iteration, true));
            t.root.joints.back().childBranches.push_back(Branch(1,t.root.joints.back().pos, iteration));
            auto &b = t.root.joints.back().childBranches.back();
            b.joints.push_back(Joint(t.pos + 3*params.ro*vec3(sin(phi),1,cos(phi)),0.1, iteration, true));
        }
    }
    else if (params.root_type == 1)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-1,0), 1, iteration, true));
        //t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0, 0.1*params.Xm*params.ro, 0), 0.9,false));
        float d = 0.05*params.Xm*params.ro;
        for (int i = 0; i < 4; i++)
        {
            t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0.1*self_rand(-1,1)*d, d, 0.1*self_rand(-1,1)*d), 0.9,iteration, false));
        }
        for (int i = 0; i < 4; i++)
        {
            t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0.1*self_rand(-1,1)*d, d, 0.1*self_rand(-1,1)*d), 0.9, iteration, true));
        }
    }
    else if (params.root_type == 2)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-2,0), 1, iteration, true));
        for (int i = 0; i < 10; i++)
        {
            float phi = 0.2*PI*(i + self_rand());
            t.root.joints.push_back(Joint(t.pos + glm::vec3(0, -2 + 0.1*(i+1)*params.ro, 0), 0.9, iteration, true));
            t.root.joints.back().childBranches.push_back(Branch(1,t.root.joints.back().pos, iteration));
            auto &b = t.root.joints.back().childBranches.back();
            b.joints.push_back(Joint(t.root.joints.back().pos + (float)self_rand(2.5,5.5)*params.ro*vec3(sin(phi),0,cos(phi)),0.1,iteration,false));
            b.joints.push_back(Joint(b.joints.back().pos + vec3(self_rand(-0.3,0.3),4,self_rand(-0.3,0.3)),0.1, iteration, true));
        }
    }
    else
    {
        logerr("GE Gen error: Unknows root type %d",params.root_type);
    }
    t.status = TreeStatus::GROWING;
}
void GETreeGenerator::set_levels_rec(Tree &t, Branch &b, GETreeParameters &params, int level)
{
    //logerr("level %d", level);
    b.level = level;
    t.max_depth = MAX(t.max_depth, level);
    for (Joint &j : b.joints)
    {
        for (Branch &br : j.childBranches)
        {
            if (br.alive)
                set_levels_rec(t, br, params, level + 1);
        }
    }
}
void GETreeGenerator::convert(Tree &src, ::Tree &dst)
{
    for (int j = 0; j < src.max_depth + 1; j++)
    {
        BranchHeap *br = new BranchHeap();
        dst.branchHeaps.push_back(br);
    }

    dst.leaves = new LeafHeap();
    dst.id = t_ids.fetch_add(1);
    dst.pos = src.pos;
    dst.type = src.type;
    dst.valid = true;

    dst.root = dst.branchHeaps[0]->new_branch();
    dst.root->type_id = dst.type->type_id;
    dst.root->self_id = ids.fetch_add(1);
    dst.root->level = 0;
    dst.root->dead = false;
    dst.root->center_self = src.pos;
    dst.root->center_par = vec3(0, 0, 0);
    dst.root->plane_coef = vec4(1, 0, 0, -src.pos.x);
    dst.root->id = dst.id;

    convert(src, dst, src.root, dst.root);

    /*
    debug("converted, branches: ");
    for (auto &bh : dst.branchHeaps)
    {
        debug("%d ",bh->branches.size());
    }
    debugnl();
    */
    
}
void GETreeGenerator::convert(Tree &src, ::Tree &dst, Branch &b_src, ::Branch *b_dst)
{
    int i = 0;
    for (auto it = b_src.joints.begin(); it != b_src.joints.end(); it++)
    {
        b_dst->joints.emplace_back();
        b_dst->joints.back().pos = it->pos;
        if (!it->leaf.edges.empty())
        {
            b_dst->joints.back().leaf = dst.leaves->new_leaf();
            b_dst->joints.back().leaf->pos = it->pos;
            b_dst->joints.back().leaf->type = 0;
            b_dst->joints.back().leaf->edges = it->leaf.edges;
        }
        if (i != 0)
        {
            auto prev = it;
            prev--;
            b_dst->segments.emplace_back();
            b_dst->segments.back().begin = prev->pos;
            b_dst->segments.back().rel_r_begin = prev->r;
            b_dst->segments.back().end = it->pos;
            b_dst->segments.back().rel_r_end = it->r;

            if (b_src.level == 2)
            {
                //b_dst->segments.back().rel_r_begin *= 4;
                //b_dst->segments.back().rel_r_end *= 4;
            }
        }
        for (Branch &chb : it->childBranches)
        {
            //if (chb.level >= 2)
            //    continue;
            //logerr("conv br %d %d",b_src.level, i);
            if (!chb.alive)
                continue;
            //logerr("conv br 3 %d %d heap %d",b_src.level, i, dst.branchHeaps.size());
            ::Branch *nb = dst.branchHeaps[chb.level]->new_branch();
            //logerr("conv br 3 %d %d",b_src.level, i);
            b_dst->joints.back().childBranches.push_back(nb);
            nb->type_id = dst.type->type_id;
            nb->self_id = ids.fetch_add(1);
            nb->level = chb.level;
            nb->dead = false;
            nb->center_self = b_dst->joints.back().pos;
            nb->center_par = b_dst->center_self;
            nb->plane_coef = vec4(1, 0, 0, -b_dst->joints.back().pos.x);
            nb->id = dst.id;
            //logerr("conv br 2 %d %d",b_src.level, i);
            convert(src, dst, chb, nb);
            if (nb->segments.front().rel_r_begin > b_dst->segments.back().rel_r_begin)
                logerr("aaaa %f %f", nb->segments.front().rel_r_begin, b_dst->segments.back().rel_r_begin);
        }
        i++;
    }
}

void GETreeGenerator::calc_light(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params)
{
    if (b.joints.empty())
        return;
    if (b.level <= params.tropism_level_base)
        b.distance_from_root = 0;
    float l = 0;
    float res = 0;
    int total_joints = b.joints.size();
    auto it = b.joints.rbegin();
    for (int i = b.joints.size() - 1; i > 0; i--)
    {
        Joint &j = *it;
        if (j.childBranches.empty())
        {
            j.light = MIN(1, SQR(1 / ( 1 + voxels.get_occlusion_trilinear(j.pos))));
            j.resource = j.light;
        }
        else
        {
            float l_lateral = 0;
            for (Branch &br : j.childBranches)
            {
                if (br.alive)
                {
                    br.distance_from_root = b.distance_from_root + i;
                    calc_light(br, voxels, params);
                    total_joints += br.total_joints;
                    l_lateral += br.total_light;
                }
            }
            float R_BH = (l + l_lateral) * (1 - params.lambda) * l_lateral /
                         (params.lambda * l + (1 - params.lambda) * l_lateral);
            j.light = l_lateral;
            j.resource = R_BH;
        }
        l += j.light;
        res += j.resource;
        it++;
    }
    b.joints.front().resource = res;
    b.joints.front().light = l;
    b.total_light = l;
    b.total_resource = res;
    b.total_joints = total_joints;
}

void GETreeGenerator::distribute_resource(Branch &b, GETreeParameters &params, float res_mult)
{
    {
        std::vector<std::pair<Joint *, float>> gp;
        auto it = b.joints.rbegin();
        for (int i = b.joints.size() - 1; i > 0; i--)
        {
            Joint &j = *it;
            bool productive = false;
            if (i == b.joints.size() - 1)
            {
                productive = true;
            }
            else if (j.childBranches.size() < params.max_branches && j.can_have_child_branches)
                productive = true;
            else
            {
                for (auto &br : j.childBranches)
                    productive = productive || b.alive;
            }
            if (productive)
            {
                gp.push_back(std::pair<Joint *, float>(&j, j.resource));
            }
            else
            {
                j.resource = 0;
            }
            it++;
        }
        int sz = gp.size();

        //insertion sort is fast for small arrays
        for (int i = 1; i < sz; i++)
            for (int j = i; j > 0 && gp[j - 1].second < gp[j].second; j--)
                std::swap(gp[j - 1], gp[j]);
        /*if (b.level <= 0)
        {
            logerr("res sum %f",b.total_resource);
            for (int i=0;i<sz;i++)
                logerr("res %f",gp[i].second);
        }*/
        
        #define W(i) (MAX(params.res_decrease_min, 1 - params.res_decrease_step * i))
        
        float denom_1 = 0, denom_2 = 0;
        float q_1 = params.res_q, q_2 = 1 - params.res_q;

        for (int i = 0; i < sz; i++)
        {
            denom_1 += gp[i].second*W(i);
            denom_2 += W(i);
        }
        denom_1 *= (q_1 + q_2);
        denom_2 *= (q_1 + q_2);
        /*
        float trm = MAX(0,params.top_res_mult_base - b.level*params.top_res_mult_level_decrease);
        b.joints.back().resource += trm*denom_1;
        denom_1 *= (1 + trm);
        if (!gp.empty())
            gp[0].second = b.joints.back().resource;
        */
        for (int i = 0; i < sz; i++)
        {
            float res = res_mult * b.total_resource * (q_1*(gp[i].second*W(i)) / denom_1 + q_2*(W(i)) / denom_2);
            auto &j = *(gp[i].first);
            int b_cnt = j.childBranches.empty() ? 1 : 0;
            for (auto &br : j.childBranches)
                b_cnt += br.alive;
            if (j.childBranches.empty())
                j.resource = res / b_cnt;
            else
                j.resource = 0;
            //if (b.level <= 0) logerr("res2 j %f",j.resource);
            for (auto &br : j.childBranches)
            {
                if (br.alive)
                {
                    br.total_resource = res / b_cnt;
                    //if (b.level <= 0) logerr("res2 br %f", br.total_resource);
                }
            }
        }
    }
    for (auto &j : b.joints)
        for (auto &br : j.childBranches)
            if (br.alive)
                distribute_resource(br, params, res_mult);
}

void GETreeGenerator::add_SPCol_points_solid_angle(vec3 pos, vec3 dir, float r_max, int cnt, float min_psi)
{
    vec3 cr, trd;
    cross_vecs(dir, cr, trd);
    float r, phi, psi;
    for (int i = 0; i < cnt; i++)
    {
        r = self_rand(0, r_max);
        phi = self_rand(0, 2 * PI);
        psi = self_rand(min_psi, PI / 2);
        vec3 dr = r * cos(psi) * sin(phi) * cr + r * cos(psi) * cos(phi) * trd + r * sin(psi) * dir;
        vec3 ps = pos + dr;
        sp_data.add(ps);
        //logerr("sp data set %f %f %f %f %f %f",dir.x, dir.y, dir.z, cr.x, cr.y, cr.z);
    }
}

void GETreeGenerator::prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params,
                                                           std::vector<GrowPoint> &growth_points,
                                                           int max_growth_per_node)
{
    float iter_frac = 1 - (float)iteration / params.max_iterations;
    int i = 0;
    for (auto it = b.joints.begin(); it != b.joints.end(); it++)
    {
        if (i == 0)
        {
            i++;
            continue;
        }
        ////logerr("%d prep",i);
        auto prev = it;
        prev--;
        auto &j = *it;
        float max_r = params.ro * max_growth_per_node;
        int sp_cnt = MAX(params.sp_points_base * iter_frac, 2);
        glm::vec3 pd = normalize(j.pos - b.joints.front().pos);
        float resource = params.resource_mult * j.resource;
        //if (b.level <= 1)
        //    resource *= 5;
        if (i == b.joints.size() - 1 && j.childBranches.size() < params.max_branches)
        {
            //grow branch forward
            GrowthType t = GrowthType::END;
            if (b.joints.size() > params.max_joints_in_branch)
                t = GrowthType::END_BRANCH;
            growth_points.push_back(GrowPoint(&j, &b, t, pd, resource,i));
            //logerr("j pos %f %f %f %f %f %f",j.pos.x,j.pos.y, j.pos.z, prev->pos.x, prev->pos.y, prev->pos.z);
            //if (sqrt(SQR(j.pos.x - g_center.x) + SQR(j.pos.z - g_center.z)) < 0.33*params.ro*params.Xm)
                add_SPCol_points_solid_angle(j.pos, pd, max_r, sp_cnt, PI / 3);
            //else if (b.level <= 2)
            //    logerr("frac %f %f",sqrt(SQR(j.pos.x - g_center.x) + SQR(j.pos.z - g_center.z)), 0.1*params.ro*params.Xm);
        }
        else if (j.can_have_child_branches && j.childBranches.size() < params.max_branches)
        {
            //create child branch
            //logerr("j pos 2 %f %f %f",j.pos.x,j.pos.y, j.pos.z);
            //if (sqrt(SQR(j.pos.x - g_center.x) + SQR(j.pos.z - g_center.z)) < 0.33*params.ro*params.Xm)
                growth_points.push_back(GrowPoint(&j, &b, GrowthType::BRANCHING, pd, resource,i));
            add_SPCol_points_solid_angle(j.pos, pd, max_r, sp_cnt, 0);
        }
        for (Branch &br : j.childBranches)
        {
            if (br.alive)
            {
                prepare_nodes_and_space_colonization(t, br, params, growth_points, max_growth_per_node);
            }
        }
        i++;
    }
}

void GETreeGenerator::grow_nodes(Tree &t, GETreeParameters &params,
                                 std::vector<GrowPoint> &growth_points,
                                 LightVoxelsCube &voxels,
                                 int max_growth_per_node)
{
    for (int i = 0; i < max_growth_per_node; i++)
    {
        std::vector<int> permutations = std::vector<int>(growth_points.size(), 0);
        for (int j = 0; j < growth_points.size(); j++)
        {
            permutations[j] = j;
        }
        std::shuffle(permutations.begin(), permutations.end(), gen);

        for (int j = 0; j < growth_points.size(); j++)
        {
            auto &gp = growth_points[permutations[j]];
            if (gp.resource_left < 1)
                gp.gType = GrowthType::FINISHED;
            //logerr("growing point %d",j);
            if (gp.gType == GrowthType::FINISHED)
                continue;
            gp.resource_left -= 1.0;

            Joint *start = gp.joint;
            Branch *br = gp.base_branch;
            glm::vec3 prev_dir = gp.prev_dir;
            float nu = 1.0;
            if (gp.gType == GrowthType::BRANCHING)
            {
                gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos, iteration));

                //choose base direction of a new branch
                vec3 b, c;
                cross_vecs(prev_dir, b, c);
                float r = params.ro;
                float phi = self_rand(-PI, PI);
                float psi = self_rand(params.branching_angle_min, params.branching_angle_max);

                vec2 planar_dir = vec2(sin(phi), cos(phi));
                planar_dir = normalize(planar_dir - 3.0f*gp.base_branch->average_chb_dir);
                gp.base_branch->average_chb_dir += planar_dir;

                prev_dir = r * cos(psi) * planar_dir.x * b + r * cos(psi) * planar_dir.y * c + r * sin(psi) * prev_dir;
                start = &(gp.joint->childBranches.back().joints.front());
                br = &(gp.joint->childBranches.back());
                br->distance_from_root = gp.base_branch->distance_from_root + gp.joint_n;
                nu *= params.branching_tropims_mult;
                prev_dir *= 3.0f;
                //logerr("new branch created");
            }
            else if (gp.gType == GrowthType::END_BRANCH)
            {
                gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos, iteration));
                start = &(gp.joint->childBranches.back().joints.front());
                br = &(gp.joint->childBranches.back());
                br->distance_from_root = gp.base_branch->distance_from_root + gp.joint_n;
            }
            float influence_r = 2 * params.ro;
            vec3 best_pos;
            float best_occ;
            if (find_best_pos(voxels, influence_r, start->pos, prev_dir, PI, best_pos, best_occ) && best_pos.x == best_pos.x)
            {
                float distance_from_root = (br->level < params.tropism_level_base) ? 0 : br->joints.size() + br->distance_from_root;
                vec3 best_dir = prev_dir + params.mu * normalize(best_pos - start->pos) +
                                (nu) * tropism(distance_from_root, params);
                vec3 new_pos = start->pos + params.ro * normalize(best_dir);
                //logerr("prev dir %f %f %f", prev_dir.x, prev_dir.y, prev_dir.z);
                //logerr("best_pos %f %f %f", best_pos.x, best_pos.y, best_pos.z);
                //logerr("joint with new pos created %f %f %f %f %f %f",new_pos.x, new_pos.y,new_pos.z,
                //                                                     best_dir.x, best_dir.y,best_dir.z);
                br->joints.push_back(Joint(new_pos, params.base_r, iteration, self_rand() < params.k));
                t.joints_total++;
                if (t.joints_total > joints_limit)
                {
                    t.status == TreeStatus::GROWN;
                    return;
                }
                //float b = params.b_min + (params.b_max - params.b_min) *
                //                             CLAMP((float)(iteration - params.tau) / (params.max_iterations - params.tau), 0, 1);
                set_occlusion_joint(br->joints.back(),1,params,voxels);
                sp_data.remove_close(new_pos, 0.75 * params.ro);

                gp.base_branch = br;
                gp.joint = &(br->joints.back());
                gp.prev_dir = normalize(br->joints.back().pos - br->joints.front().pos);
                gp.gType = GrowthType::END;
            }
            else
            {
                //logerr("failed to find best pos");
                //we cannot grow this branch
                gp.gType = GrowthType::FINISHED;
            }
        }
    }
}

void GETreeGenerator::remove_branches(Tree &t, Branch &b, GETreeParameters &params, LightVoxelsCube &voxels)
{
    b.total_resource = 0;
    if (b.level < 2)
    {
        //logerr("level %d joints %d params %f %f %f",b.level, b.joints.size(), params.r_s, params.branching_angle_max, 
        //       params.branching_angle_min);
    }
    if (b.joints.size() < 2)
    {
        b.total_joints = b.joints.size();
        b.total_resource = 0;
        return;
    }
    float total_light = 0;
    int total_joints = b.joints.size();
    int i = 0;
    float dead_b_count = 0;
    float ch_b_count = 0;
    for (Joint &j : b.joints)
    {
        if (i > 0)
        {
            if (j.childBranches.empty())
                total_light += MIN(1, 1 / (0.5 + MAX(voxels.get_occlusion_trilinear(j.pos),0)));
            if (true)
            {
                //logerr("%f %f aaa", voxels.get_occlusion_trilinear(j.pos),  MIN(1, 1 / (0.5 + MAX(voxels.get_occlusion_trilinear(j.pos),0))));
            }
            //logerr("l %f", voxels.get_occlusion_trilinear(j.pos));
            for (Branch &br : j.childBranches)
            {
                if (br.alive)
                {
                    remove_branches(t, br, params, voxels);
                    ch_b_count++;
                    total_light += br.total_resource; //total_light
                    if (br.total_resource != br.total_resource)
                    {
                        logerr("nan detected here");
                    }
                    total_joints += br.total_joints;
                }
                else
                    dead_b_count++;
            }

            if (true)
            {
                auto it = j.childBranches.begin();
                while (it != j.childBranches.end())
                {
                    if (it->alive || it->base_r >= 3*params.base_r)
                    {
                        it++;
                    }
                    else
                    {
                        it = j.childBranches.erase(it);
                    }
                }
            }
        }
        i++;
    }
    float lq = total_light / total_joints + 1e-6;
    float size_log = MIN(log2f(total_joints), 10);
    float remove_q = lq/(params.r_s + params.rs_size_factor * size_log);
    float q = remove_q > 1 ? 0.5*size_log*(remove_q - 1) : 0.5*size_log*(1/remove_q - 1);
    float remove_chance = 1 - 1/(1 + exp(-q)); 
    //if (params.root_type == 1)
    //    logerr("params.r_s %d %f %f chance = %f %f %f", total_joints, total_light, params.r_s + params.rs_size_factor * MIN(log2f(total_joints), 10),
    //          remove_q, q, remove_chance);
    //logerr("%f tot", total_light);
    if (total_light != total_light)
    {
        logerr("NAN detected");
        int r =  0;
        for (Joint &j : b.joints)
        {   
            for (Branch &br : j.childBranches)
            {
                logerr("br %d res %f", r, br.total_resource);
                r++;
            }
        }
    }
    if (total_joints > 10)
    {
        //logerr("br is ok %d %d %d %f",b.level, total_joints, b.joints.size(), ch_b_count);
    }
    if (b.level >= params.remove_min_level && (total_joints > 10) && 
       (self_rand() < remove_chance ||
        total_joints > 50 && (float)b.joints.size()/total_joints > 0.33 ||
        dead_b_count/b.joints.size() > 0.75 ||
        b.joints.size() > 0.5*params.max_joints_in_branch && ch_b_count < 2))
    {
        //remove branch
        //logerr("remove branch %d %f",total_joints, total_light);
        set_occlusion(b, voxels, params, -1); //remove occlusion from deleted branch
        b.joints = {};
        b.alive = false;
        b.total_joints = 0;
        b.total_resource = 0;
    }
    else
    {
        b.total_resource = total_light;
        b.total_joints = total_joints;
    }
}

void GETreeGenerator::recalculate_radii(Tree &t, Branch &b, GETreeParameters &params)
{
    float r = pow(params.base_r, params.r_pow);
    for (auto it = b.joints.rbegin(); it != b.joints.rend(); it++)
    {
        for (Branch &br : it->childBranches)
        {
            if (br.alive)
            {
                recalculate_radii(t, br, params);
            }
            else
            {
                br.base_r *= 0.75;
            }
            r += pow(br.base_r, params.r_pow); //use memory of removed branches too
        }
        it->r = pow(r, 1.0 / params.r_pow);
    }
    b.base_r = pow(r, 1.0 / params.r_pow);
}

void GETreeGenerator::set_occlusion(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params, float mul)
{
    int i = 0;
    for (Joint &j : b.joints)
    {
        if (i > 0)
        {
            set_occlusion_joint(j, mul, params, voxels);
            for (Branch &br : j.childBranches)
            {
                if (br.alive)
                    set_occlusion(br, voxels, params, mul);
            }
        }
        i++;
    }
}

void GETreeGenerator::set_occlusion_joint(Joint &j, float base_value, GETreeParameters &params, LightVoxelsCube &voxels)
{
    glm::vec3 p = j.pos;
    int rnd_seed = 17*abs(p.x) + 19*abs(p.y) + 23*abs(p.z);
    voxels.set_occluder_pyramid_fast(j.pos, base_value, params.occlusion_pyramid_d, rnd_seed);
}

bool GETreeGenerator::find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos, glm::vec3 dir, float angle,
                                    glm::vec3 &best_pos, float &best_occ)
{
    return sp_data.find_best_pos(voxels, r, pos, dir, angle, best_pos, best_occ);
}

bool GETreeGenerator::SpaceColonizationData::find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos,
                                                           glm::vec3 dir, float angle,
                                                           glm::vec3 &best_pos, float &best_occ)
{
    if (!active)
        return false;
    best_occ = 1000;
    float cs = cos(angle);

    std::function<void(glm::vec3 &)> func = [&](glm::vec3 &p)
    {
        if (dot(normalize(p - pos), dir) > cs)
        {
            float occ = voxels.get_occlusion_trilinear(p);
            //logerr("%f occ", occ);
            if (occ < best_occ)
            {
                best_occ = occ;
                best_pos = p;
            }
        }
    };

    AABB box = AABB(pos - r * vec3(1, 1, 1), pos + r * vec3(1, 1, 1));
    octree.apply_to_neighbours_sphere(box, r, pos, func);

    if (best_occ >= 1000)
        return false;
    else
        return true;
}

void GETreeGenerator::SpaceColonizationData::remove_close(glm::vec3 pos, float r)
{
    if (!active)
       return;
    AABB box = AABB(pos - r * vec3(1, 1, 1), pos + r * vec3(1, 1, 1));
    octree.remove_in_sphere(box, r, pos);
}

void GETreeGenerator::SpaceColonizationData::add(glm::vec3 pos)
{
    if (!active)
        return;
    positions.push_back(pos);
}

void GETreeGenerator::SpaceColonizationData::prepare(LightVoxelsCube &voxels)
{
    if (!active)
        return;
    octree.create(voxels.get_bbox());
    octree.insert_vector(positions);
    positions.clear();
}
