#include "GE_generator.h"
#include <algorithm>
using namespace glm;

GETreeParameters GETreeGenerator::defaultParameters = GETreeParameters();

//TreeTypeData def_ttd = TreeTypeData(-1,&(GETreeGenerator::defaultParameters), "wood","leaf");
bool GETreeGenerator::iterate(LightVoxelsCube &voxels)
{
    bool growing = false;
    for (auto &t : trees)
    {
        const GETreeParameters *p_ptr = dynamic_cast<const GETreeParameters *>(t.type->get_params());
        const GETreeParameters &params = p_ptr ? *(p_ptr) : defaultParameters;
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
            int max_growth = MIN(32, round(dX_dt));

            //logerr("iteration %d max growth %d", iteration, max_growth);
            if (iteration >= params.max_iterations || max_growth <= 0)
                t.status = TreeStatus::GROWN;
            else
            {
                std::vector<GrowPoint> growth_points;
                sp_data = SpaceColonizationData();
                sp_data.active = true;
                calc_light(t.root, voxels, params);
                float l0 = t.root.total_light;
                t.root.total_light *= MAX((t.root.total_light + 30)/(t.root.total_light + 0.1),1.5);
                float l1 = CLAMP(t.root.total_light, 0, 15*params.Xm);
                float l2 = CLAMP(t.root.total_light - l1, 0, 30*params.Xm);
                float l3 = CLAMP(t.root.total_light - l2, 0, 60*params.Xm);
                float l4 = CLAMP(t.root.total_light - l2, 0, 1000*params.Xm);
                float l_corr = l1 + 0.4*l2 + 0.16*l3 + 0.0*l4;
                
                calc_distance_from_root(t.root, params, glm::vec2(0,0));
                distribute_resource(t.root, params, 1.0);
                prepare_nodes_and_space_colonization(t, t.root, params, growth_points, max_growth);
                sp_data.prepare(voxels);
                grow_nodes(t, params, growth_points, voxels, max_growth);
                remove_branches(t, t.root, params, voxels);

                t.iteration++;
                growing = true;
            }
        }
    }

    return growing;
}
void GETreeGenerator::plant_tree(glm::vec3 pos, const TreeTypeData *type)
{
    GETreeParameters *ps = dynamic_cast<GETreeParameters *>(type->get_params());
    if (!ps)
    {
        logerr("Tree type %d cannot be generated with GE tree generator", type->type_id);
        return;
        //type->params = &(GETreeGenerator::defaultParameters);
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
        GETreeParameters *p_ptr = dynamic_cast<GETreeParameters *>(t.type->get_params());
        const GETreeParameters &params = p_ptr ? *(p_ptr) : defaultParameters;
        set_levels_rec(t, t.root, params, 0);
        recalculate_radii(t, t.root, params);
        create_leaves(t.root, params, 2, voxels);
        convert(t, trees_external[trees.size() - i - 1]);
        i++;
    }
    trees.clear();
    iteration = 0;
}

void GETreeGenerator::create_leaves(Branch &b, const GETreeParameters &params, int level_from, LightVoxelsCube &voxels)
{
    if (!b.alive || b.joints.size() < 2)
        return;
    auto jit = b.joints.begin();
    jit++;
    auto prev_jit = b.joints.begin();
    while (jit != b.joints.end())
    {
        Joint &j = *jit;
        vec3 dir_x = normalize(jit->pos - prev_jit->pos);
        vec3 dir_y, dir_z;
        cross_vecs(dir_x, dir_y, dir_z);

        if (j.childBranches.empty() && b.level >= level_from && params.leaves_cnt > 0 &&
            j.r < params.base_r * params.leaves_max_r)
        {
            float f_leaves_cnt = params.leaves_cnt * SQR(1 / (0.5 + voxels.get_occlusion_trilinear(j.pos)));
            int l_cnt = f_leaves_cnt;
            if (self_rand() < (f_leaves_cnt - l_cnt))
                l_cnt++;
            
            for (int i=0;i<l_cnt;i++)
            {
                glm::vec3 rd1 = normalize(vec3(self_rand(-params.leaves_angle_a, params.leaves_angle_a), 
                                               self_rand(-1, 1), self_rand(-1, 1)));
                rd1 = dir_x*rd1.x + dir_y*rd1.y + dir_z*rd1.z;

                glm::vec3 rd2 = normalize(vec3(self_rand(-params.leaves_angle_b, params.leaves_angle_b),
                                               self_rand(-1, 1), self_rand(-1, 1)));
                rd2 = dir_x*rd2.x + dir_y*rd2.y + dir_z*rd2.z;
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
        jit++;
        prev_jit++;
    }
}
void GETreeGenerator::create_initial_trunk(Tree &t, const GETreeParameters &params)
{
    long param_dependant_seed = 71*params.initial_trunk_scale.x  + 
                                73*params.initial_trunk_scale.y  +
                                79*abs(params.trunk_bonus_radius)+
                                173*params.trunk_bonus_radius_mod;
    for (int i=0;i<param_dependant_seed % 10; i++)
        self_rand();                            
    t.root = Branch();
    t.root.level = 0;
    t.root.can_be_removed = false;
    float sz_x = params.initial_trunk_scale.x*params.ro;
    float sz_y = params.initial_trunk_scale.y*params.ro;
    if (params.root_type == 0)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-1,0), iteration, true));
        for (int i = 0; i < 10; i++)
        {
            float phi = 0.2*PI*(i + self_rand());
            t.root.joints.push_back(Joint(t.pos + glm::vec3(0, 0.1*(i+1), 0), iteration, true));
            t.root.joints.back().childBranches.push_back(Branch(1,t.root.joints.back().pos, iteration,false));
            auto &b = t.root.joints.back().childBranches.back();
            b.joints.push_back(Joint(t.pos + 3.0f*vec3(sz_x*sin(phi),sz_y,sz_x*cos(phi)), iteration, true));
        }
    }
    else if (params.root_type == 1)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-1,0), iteration, true, params.trunk_bonus_radius));
        float dx = 0.05*params.Xm*sz_x;
        float dy = 0.05*params.Xm*sz_y;
        for (int i = 0; i < 4; i++)
        {
            t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0.1*self_rand(-1,1)*dx, dy, 0.1*self_rand(-1,1)*dx), 
                                          iteration, false, params.trunk_bonus_radius*pow(1 - (i+1)/9.0, params.trunk_bonus_radius_mod)));
        }
        for (int i = 0; i < 4; i++)
        {
            t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0.1*self_rand(-1,1)*dx, dy, 0.1*self_rand(-1,1)*dx), 
                                          iteration, true, params.trunk_bonus_radius*pow(1 - (i+5)/9.0, params.trunk_bonus_radius_mod)));
        }
    }
    else if (params.root_type == 2)
    {
        t.root.joints.push_back(Joint(t.pos + glm::vec3(0,-2,0), iteration, true));
        for (int i = 0; i < 10; i++)
        {
            float phi = 0.2*PI*(i + self_rand());
            t.root.joints.push_back(Joint(t.pos + glm::vec3(0, (-2 + 0.03*(i+1))*sz_y, 0), iteration, true));
            t.root.joints.back().childBranches.push_back(Branch(1,t.root.joints.back().pos, iteration,false));
            auto &b = t.root.joints.back().childBranches.back();
            b.joints.push_back(Joint(t.root.joints.back().pos + (float)self_rand(1,5)*sz_x*vec3(sin(phi),0,cos(phi)), iteration,false));
            b.joints.push_back(Joint(b.joints.back().pos + vec3(self_rand(-0.2,0.2),4*sz_y,self_rand(-0.2,0.2)), iteration, true));
        }
    }
    else if (params.root_type == 3)
        create_universal_initial_trink(t, params);
    else
    {
        logerr("GE Gen error: Unknows root type %d",params.root_type);
    }
    t.status = TreeStatus::GROWING;
}

void GETreeGenerator::create_universal_initial_trink(Tree &t, const GETreeParameters &params)
{
    t.status = TreeStatus::GROWING;
    t.root = Branch();
    t.root.level = 0;
    t.root.can_be_removed = false;
    float sz_x = params.initial_trunk_scale.x*params.ro;
    float sz_y = params.initial_trunk_scale.y*params.ro;

    float dL;
    if (params.L0 > 0)
    {
        t.root.joints.push_back(Joint(t.pos, 0, false, params.R0));
        dL = sz_y*params.L0/4;
    }
    else
    {
        t.root.joints.push_back(Joint(t.pos - glm::vec3(0,sz_y*params.L0,0), 0, false, params.R0));
        dL = -sz_y*params.L0/4;        
    }

    for (int i=0;i<4;i++)
    {
        glm::vec3 pos = t.root.joints.front().pos + glm::vec3(0,(i+1)*dL,0);
        pos.x += params.dX*sz_x*self_rand(-1,1);
        pos.z += params.dX*sz_x*self_rand(-1,1);
        float r = params.iR[i] * t.root.joints.back().r;
        t.root.joints.push_back(Joint(pos, 0, false, r));
    }

    int prev_pos_phi_n = 0;
    for (int i=0;i<params.B_cnt;i++)
    {
        std::vector<int> pos_phi_n;
        if (i < 4)
            pos_phi_n.push_back(i);
        else 
        {
            for (int j=0;j<4;j++)
            {
                if (j != prev_pos_phi_n)
                    pos_phi_n.push_back(j);
            }
        }
        #define RND_N (i < 4 ? i : ((int)self_rand(0,3.999)))

        glm::vec3 r_pos = t.root.joints.back().pos;
        r_pos.y += sz_y*params.Lt[RND_N];
        r_pos.x = t.pos.x + params.dX*sz_x*self_rand(-1,1);
        r_pos.z = t.pos.z + params.dX*sz_x*self_rand(-1,1);

        float r = t.root.joints.back().r * params.dR[RND_N];

        float sectors = 8;
        float phi = ((int)(sectors*params.phi[pos_phi_n[(int)self_rand(0,pos_phi_n.size() - 1e-4)]]) + self_rand(0,1))/sectors*PI;
        float psi = ((int)(sectors*params.psi[RND_N]) + self_rand(0,1))/sectors*PI/2;
        float l = params.Lb[RND_N];

        glm::vec3 b_pos = r_pos + l*glm::vec3(sz_x*sin(psi)*cos(phi),sz_y*cos(psi),sz_x*sin(psi)*sin(phi));

        t.root.joints.push_back(Joint(r_pos,0, true, r));
        t.root.joints.back().childBranches.push_back(Branch(1, r_pos, 0, false));
        auto &b = t.root.joints.back().childBranches.back();
        b.joints.push_back(Joint(b_pos, 0, true, 0));
    }
}

void GETreeGenerator::set_levels_rec(Tree &t, Branch &b, const GETreeParameters &params, int level)
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
    dst.id = tree_next_id.fetch_add(1);
    dst.pos = src.pos;
    dst.type = src.type;
    dst.valid = true;

    dst.root = dst.branchHeaps[0]->new_branch();
    dst.root->type_id = dst.type->type_id;
    dst.root->self_id = branch_next_id.fetch_add(1);
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
            nb->self_id = branch_next_id.fetch_add(1);
            nb->level = chb.level;
            nb->dead = false;
            nb->center_self = b_dst->joints.back().pos;
            nb->center_par = b_dst->center_self;
            nb->plane_coef = vec4(1, 0, 0, -b_dst->joints.back().pos.x);
            nb->id = dst.id;
            //logerr("conv br 2 %d %d",b_src.level, i);
            convert(src, dst, chb, nb);
            //if (nb->segments.front().rel_r_begin > b_dst->segments.back().rel_r_begin)
            //    logerr("aaaa %f %f", nb->segments.front().rel_r_begin, b_dst->segments.back().rel_r_begin);
        }
        i++;
    }
}

void GETreeGenerator::calc_distance_from_root(Branch &b, const GETreeParameters &params, glm::vec2 start_distance)
{
    if (b.joints.size() < 2)
        return;
    auto it = b.joints.begin();
    it->distance_from_root = start_distance;
    it++;
    auto prev_it = b.joints.end();
    while (it != b.joints.end())
    {
        vec3 d = it->pos - prev_it->pos;
        it->distance_from_root = prev_it->distance_from_root + vec2(sqrt(SQR(d.x) + SQR(d.z)), abs(d.y));
        for (auto &ch_b : it->childBranches)
        {
            calc_distance_from_root(ch_b, params, it->distance_from_root);
        }
        it++;
        prev_it++;
    }
}

void GETreeGenerator::calc_light(Branch &b, LightVoxelsCube &voxels, const GETreeParameters &params)
{
    if (b.joints.empty())
        return;
    if (b.level <= params.tropism_level_base)
        b.distance_from_trunk = 0;
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
                    br.distance_from_trunk = b.distance_from_trunk + i;
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

void GETreeGenerator::distribute_resource(Branch &b, const GETreeParameters &params, float res_mult)
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

        //bubble sort is fast for small arrays
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
        #define RD(d) (CLAMP(1 - (d.x-params.distance_fine_start.x)/(params.distance_fine_start.x + params.distance_fine_slope.x),0,1)* \
                       CLAMP(1 - (d.y-params.distance_fine_start.y)/(params.distance_fine_start.y + params.distance_fine_slope.y),0,1))
        float denom_1 = 0, denom_2 = 0;
        float q_1 = params.res_q, q_2 = 1 - params.res_q;

        for (int i = 0; i < sz; i++)
        {
            denom_1 += gp[i].second*W(i)*RD(gp[i].first->distance_from_root);
            denom_2 += W(i)*RD(gp[i].first->distance_from_root);
        }

        for (int i = 0; i < sz; i++)
        {
            float res = res_mult * b.total_resource * (q_1*(gp[i].second*W(i)*RD(gp[i].first->distance_from_root)) / denom_1 + 
                                                       q_2*(W(i)*RD(gp[i].first->distance_from_root)) / denom_2);
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

void GETreeGenerator::prepare_nodes_and_space_colonization(Tree &t, Branch &b, const GETreeParameters &params,
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

void GETreeGenerator::grow_nodes(Tree &t, const GETreeParameters &params,
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
                gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos, iteration, true));

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
                br->distance_from_trunk = gp.base_branch->distance_from_trunk + gp.joint_n;
                nu *= params.branching_tropims_mult;
                prev_dir *= 3.0f;
                //logerr("new branch created");
            }
            else if (gp.gType == GrowthType::END_BRANCH)
            {
                gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos, iteration, true));
                start = &(gp.joint->childBranches.back().joints.front());
                br = &(gp.joint->childBranches.back());
                br->distance_from_trunk = gp.base_branch->distance_from_trunk + gp.joint_n;
            }
            float influence_r = 2 * params.ro;
            vec3 best_pos;
            float best_occ;
            if (find_best_pos(voxels, influence_r, start->pos, prev_dir, PI, best_pos, best_occ) && best_pos.x == best_pos.x)
            {
                float distance_from_trunk = (br->level < params.tropism_level_base) ? 0 : br->joints.size() + br->distance_from_trunk;
                vec3 best_dir = prev_dir + params.mu * normalize(best_pos - start->pos) +
                                (nu) * tropism(distance_from_trunk, params);
                vec3 new_pos = start->pos + params.ro * normalize(best_dir);
                //logerr("prev dir %f %f %f", prev_dir.x, prev_dir.y, prev_dir.z);
                //logerr("best_pos %f %f %f", best_pos.x, best_pos.y, best_pos.z);
                //logerr("joint with new pos created %f %f %f %f %f %f",new_pos.x, new_pos.y,new_pos.z,
                //                                                     best_dir.x, best_dir.y,best_dir.z);
                br->joints.push_back(Joint(new_pos, iteration, self_rand() < params.k));
                t.joints_total++;
                if (t.joints_total > joints_limit)
                {
                    t.status = TreeStatus::GROWN;
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

void GETreeGenerator::remove_branches(Tree &t, Branch &b, const GETreeParameters &params, LightVoxelsCube &voxels)
{
    b.total_resource = 0;

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
                total_light += MIN(1, 1 / (0.5 + MAX(voxels.get_occlusion_simple(j.pos),0)));
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

            auto it = j.childBranches.begin();
            while (it != j.childBranches.end())
            {
                if (it->alive)
                    it++;
                else
                    it = j.childBranches.erase(it);
            }
        }
        i++;
    }
    float lq = total_light / total_joints + 1e-6;
    float size_log = MIN(log2f(total_joints), 10);
    float remove_q = lq/(params.r_s + params.rs_size_factor * size_log);
    float q = remove_q > 1 ? 0.5*size_log*(remove_q - 1) : 0.5*size_log*(1/remove_q - 1);
    float remove_chance = 1 - 1/(1 + exp(-q)); 

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

    if (b.level >= params.remove_min_level && (total_joints > 10) && b.can_be_removed &&
       (self_rand() < remove_chance ||
        (total_joints > 50 && (float)b.joints.size()/total_joints > 0.33) ||
        dead_b_count/b.joints.size() > 0.75 ||
        (b.joints.size() > 0.5*params.max_joints_in_branch && ch_b_count < 2)))
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

void GETreeGenerator::recalculate_radii(Tree &t, Branch &b, const GETreeParameters &params)
{
    //it should be called only once - during finalize step
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
        //previously it->r contains bonus r from generation step
        it->r = CLAMP(it->r + pow(r, 1.0 / params.r_pow), params.base_r, 1000*params.base_r);
    }
    b.base_r = b.joints.front().r;
}

void GETreeGenerator::set_occlusion(Branch &b, LightVoxelsCube &voxels, const GETreeParameters &params, float mul)
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

void GETreeGenerator::set_occlusion_joint(Joint &j, float base_value, const GETreeParameters &params, LightVoxelsCube &voxels)
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
