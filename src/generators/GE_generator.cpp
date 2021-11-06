#include "GE_generator.h"
#include <algorithm>
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
    ivec3 voxels_sizes = ivec3(2*params.Xm, 3.5*params.Xm, 2*params.Xm);
    LightVoxelsCube voxels = LightVoxelsCube(t.pos, voxels_sizes, 0.5*params.ro);
    set_occlusion(t.root, voxels, params, 1);

    for (int i=0;i<params.max_iterations;i++)
    {
        iteration = i;

        auto &p = params;
        float A = ((float)p.Xm/p.X0 - 1)*exp(-p.r*iteration);
        float dX_dt = p.Xm*p.r*A/SQR(1 + A);
        int max_growth = round(dX_dt);

        logerr("iteration %d max growth %d",i, max_growth);        
        if (max_growth <= 0)
            break;
        //max_growth = 1;
        std::vector<GrowPoint> growth_points;
        SpaceColonizationData sp_data;

        calc_light(t.root, voxels, params);
        distribute_resource(t.root, params);
        prepare_nodes_and_space_colonization(t, t.root, params, growth_points, sp_data, max_growth);
        sp_data.prepare(voxels);
        //logerr("prepared %d grow points %d sp dots", growth_points.size(), sp_data.positions.size());
        grow_nodes(t, params, growth_points, sp_data, voxels, max_growth);
        recalculate_radii(t, t.root, params);
        remove_branches(t, t.root, params, voxels);
    }
    //generate tree
    //LightVoxelsCube voxels = LightVoxelsCube(t.pos, voxels_sizes, 0.5*params.ro);
    //set_occlusion(t.root, voxels, params);
    set_levels_rec(t, t.root,params,0);
    create_leaves(t.root, params, 0, voxels);
}
void GETreeGenerator::create_leaves(Branch &b, GETreeParameters &params, int level_from, LightVoxelsCube &voxels)
{
    for (Joint &j : b.joints)
    {
        if (j.childBranches.empty() && b.level >= level_from && params.leaves_cnt > 0 && 
            j.r < params.base_r*params.leaves_max_r &&
            urand() < params.leaves_cnt * (1/(0.5+voxels.get_occlusion(j.pos))))
        {
            glm::vec3 rd1 = normalize(vec3(urand(-1,1),urand(-1,1),urand(-1,1)));
            glm::vec3 rd2 = normalize(vec3(urand(-1,1),urand(-1,1),urand(-1,1)));
            float sz = params.ro*params.leaf_size_mult;
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
        for (Branch &br : j.childBranches)
        {
            if (br.alive)
                create_leaves(br, params, level_from, voxels);
        }
    }
}
void GETreeGenerator::create_initial_trunk(Tree &t, GETreeParameters &params)
{
    t.root.level = 0;
    t.root.joints.push_back(Joint(t.pos, 1));
    t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0,10*params.ro,0), 0.9));
    for (int i=0;i<10;i++)
        t.root.joints.push_back(Joint(t.root.joints.back().pos + glm::vec3(0,3*params.ro,0), 0.9));

    return;

    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,60,0), 0.7));

    t.root.joints.back().childBranches.emplace_back();
    auto &b = t.root.joints.back().childBranches.back();
    vec3 p = t.root.joints.back().pos;
    b.joints.push_back(Joint(p,0.3));
    b.joints.push_back(Joint(p + vec3(10,10,0),0.3));
    b.joints.push_back(Joint(p + vec3(13,17,0),0.3));

    b.joints.back().childBranches.emplace_back();
    auto &b1 = b.joints.back().childBranches.back();
    p = b.joints.back().pos;

    b.joints.push_back(Joint(p + vec3(16,25,0),0.3));

    b1.joints.push_back(Joint(p,0.1));
    b1.joints.push_back(Joint(p + vec3(0,05,05),0.1));
    b1.joints.push_back(Joint(p + vec3(0,10,06),0.1));
    b1.joints.push_back(Joint(p + vec3(03,11,07),0.1));

    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,90,0), 0.53));
    t.root.joints.push_back(Joint(t.pos + glm::vec3(0,120,0), 0.5));
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

            if (b_src.level == 0)
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
            nb->type_id = 0;
            nb->self_id = ids++;
            nb->level = chb.level;
            nb->dead = false;
            nb->center_self = b_dst->joints.back().pos;
            nb->center_par = b_dst->center_self;
            nb->plane_coef = vec4(1,0,0,-b_dst->joints.back().pos.x);
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
        float l = 0;
        float res = 0;
        int total_joints = b.joints.size();
        auto it = b.joints.rbegin();
        for (int i=b.joints.size()-1;i>0;i--)
        {
            Joint &j = *it;
            if (j.childBranches.empty())
            {
                j.light = MIN(1, 1/(0.5 + voxels.get_occlusion_simple(j.pos)));
                j.resource = j.light;
            }
            else
            {
                float l_lateral = 0;
                for (Branch &br : j.childBranches)
                {
                    if (br.alive)
                    {
                        calc_light(br, voxels, params);
                        total_joints += br.total_joints;
                        l_lateral += br.total_light;
                    }
                }
                float R_BH = (l + l_lateral)*(1-params.lambda)*l_lateral/
                             (params.lambda*l + (1-params.lambda)*l_lateral);
                j.light = l_lateral;
                j.resource = R_BH;
            }
            l+= j.light;
            res+= j.resource;
            it++;
        }
        b.joints.front().resource = res;
        b.joints.front().light = l;
        b.total_light = l;
        b.total_resource = res;
        b.total_joints = total_joints;
    }

    void GETreeGenerator::distribute_resource(Branch &b, GETreeParameters &params)
    {
        {
        std::vector<std::pair<Joint *,float>> gp;
        auto it = b.joints.rbegin();
        for (int i=b.joints.size()-1;i>0;i--)
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
                gp.push_back(std::pair<Joint *,float>(&j,j.resource));
            }
            else
            {
                j.resource = 0;
            }
            it++;
        }
        int sz = gp.size();
        
        //insertion sort is fast for small arrays
        for(int i=1;i<sz;i++)     
	        for(int j=i;j>0 && gp[j-1].second<gp[j].second;j--)
			    std::swap(gp[j-1],gp[j]);
       /*if (b.level <= 0)
        {
            logerr("res sum %f",b.total_resource);
            for (int i=0;i<sz;i++)
                logerr("res %f",gp[i].second);
        }*/
        float denom = 0;
        #define W(i) MAX(0.03,1-0.07*i)
        for (int i=0;i<sz;i++)
            denom += W(i);
        /*float trm = MAX(0,params.top_res_mult_base - b.level*params.top_res_mult_level_decrease);
        b.joints.back().resource += trm*denom;
        denom *= (1 + trm);*/
        for (int i=0;i<sz;i++)
        {
            float res = b.total_resource*(W(i))/denom;
            auto &j = *(gp[i].first);
            int b_cnt = j.childBranches.empty() ? 1 : 0;
            for (auto &br : j.childBranches)
                b_cnt += br.alive;
            if (j.childBranches.empty())
                j.resource = res/b_cnt;
            else
                j.resource = 0;
            //if (b.level <= 0) logerr("res2 j %f",j.resource);
            for (auto &br : j.childBranches)
            {
                if (br.alive)
                {
                    br.total_resource = res/b_cnt;
                    //if (b.level <= 0) logerr("res2 br %f", br.total_resource);
                }
                
            }
        }
        }
        for (auto &j : b.joints)
            for (auto &br : j.childBranches)
                if (br.alive)
                    distribute_resource(br, params);
    }

    void cross_vecs(vec3 a, vec3 &b, vec3 &c)
    {
        b = vec3(1,0,0);
        if (abs(dot(b-a,b-a)) > 1 - 1e-6)
            b = vec3(0,0,1);
        b = cross(a, b);
        c = cross(a, b);
    }
    void GETreeGenerator::add_SPCol_points_solid_angle(vec3 pos, vec3 dir, float r_max, int cnt, float min_psi, 
                                                       SpaceColonizationData &sp_data)
    {
        vec3 cr, trd;
        cross_vecs(dir, cr, trd);
        float r, phi, psi;
        for (int i=0;i<cnt;i++)
        {
            r = urand(0, r_max);
            phi = urand(0, 2*PI);
            psi = urand(min_psi, PI/2);
            vec3 dr = r*cos(psi)*sin(phi)*cr + r*cos(psi)*cos(phi)*trd + r*sin(psi)*dir;
            vec3 ps = pos + dr;
            sp_data.add(ps);
            //logerr("sp data set %f %f %f %f %f %f",dir.x, dir.y, dir.z, cr.x, cr.y, cr.z);
        }
    }

    void GETreeGenerator::prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params, 
                                              std::vector<GrowPoint> &growth_points,
                                              SpaceColonizationData &sp_data,
                                              int max_growth_per_node)
    {
        float iter_frac = 1 - (float)iteration/params.max_iterations;
        int i=0;
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
            float max_r = 2*params.ro * max_growth_per_node;
            int sp_cnt = MAX(params.sp_points_base * iter_frac + 15, 2);
            glm::vec3 pd = normalize(j.pos - b.joints.front().pos);
            float resource = params.r*params.resource_mult*j.resource;
            //if (b.level <= 1)
            //    resource *= 5;
            if (i == b.joints.size() - 1)
            {
                //grow branch forward
                GrowthType t = GrowthType::END;
                //if (b.joints.size() > params.max_joints_in_branch)
                //    t = GrowthType::END_BRANCH;
                growth_points.push_back(GrowPoint(&j, &b, t, pd, resource));
                //logerr("j pos %f %f %f %f %f %f",j.pos.x,j.pos.y, j.pos.z, prev->pos.x, prev->pos.y, prev->pos.z);
                add_SPCol_points_solid_angle(j.pos, pd, max_r,sp_cnt,PI/3,sp_data);
            }
            else if (j.can_have_child_branches && j.childBranches.size() < params.max_branches)
            {
                //create child branch
                //logerr("j pos 2 %f %f %f",j.pos.x,j.pos.y, j.pos.z);
                growth_points.push_back(GrowPoint(&j, &b, GrowthType::BRANCHING, pd, resource));
                add_SPCol_points_solid_angle(j.pos, pd, max_r,sp_cnt,0,sp_data);
            }
            for (Branch &br : j.childBranches)
            {
                if (br.alive)
                {
                    prepare_nodes_and_space_colonization(t,br,params,growth_points,sp_data,max_growth_per_node);
                }
            }
            i++;
        }
    }

    vec3 tropism(float n, GETreeParameters &params)
    {
        return vec3(0,MAX(0.5-SQR(3*n/params.Xm),-1),0);
    }
    
    void GETreeGenerator::grow_nodes(Tree &t, GETreeParameters &params, 
                    std::vector<GrowPoint> &growth_points,
                    SpaceColonizationData &sp_data,
                    LightVoxelsCube &voxels,
                    int max_growth_per_node)
    {
        for (int i=0;i<max_growth_per_node;i++)
        {
            std::vector<int> permutations = std::vector<int>(growth_points.size(),0);
            for (int j=0;j<growth_points.size();j++)
            {
                permutations[j] = j;
            }
            std::random_shuffle(permutations.begin(), permutations.end());


            for (int j=0;j<growth_points.size();j++)
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
                if (gp.gType == GrowthType::BRANCHING)
                {
                    gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos));

                    //choose base direction of a new branch
                    vec3 b,c;
                    cross_vecs(prev_dir, b, c);
                    float r = params.ro;
                    float phi = urand(0, 2*PI);
                    float psi = urand(params.branching_angle_min,params.branching_angle_max);
                    
                    prev_dir = r*cos(psi)*sin(phi)*b + r*cos(psi)*cos(phi)*c + r*sin(psi)*prev_dir;
                    start = &(gp.joint->childBranches.back().joints.front()); 
                    br = &(gp.joint->childBranches.back());
                    //logerr("new branch created");
                }
                else if (gp.gType == GrowthType::END_BRANCH)
                {
                    gp.joint->childBranches.push_back(Branch(gp.base_branch->level + 1, gp.joint->pos));
                    start = &(gp.joint->childBranches.back().joints.front()); 
                    br = &(gp.joint->childBranches.back());
                }
                float influence_r = 4*params.ro;
                vec3 best_pos;
                float best_occ;
                if (sp_data.find_best_pos(voxels, influence_r, start->pos, prev_dir, PI/6, best_pos, best_occ) 
                    && best_pos.x == best_pos.x)
                {
                    vec3 best_dir = prev_dir + params.mu*normalize(best_pos - start->pos) + 
                                    (params.nu)*tropism(br->joints.size(),params);
                    vec3 new_pos =  start->pos +  params.ro*normalize(best_dir);
                    //logerr("prev dir %f %f %f", prev_dir.x, prev_dir.y, prev_dir.z);
                    //logerr("best_pos %f %f %f", best_pos.x, best_pos.y, best_pos.z);
                    //logerr("joint with new pos created %f %f %f %f %f %f",new_pos.x, new_pos.y,new_pos.z,
                    //                                                     best_dir.x, best_dir.y,best_dir.z);
                    br->joints.push_back(Joint(new_pos, params.base_r, urand() < params.k));

                    float b = params.b_min + (params.b_max - params.b_min)*
                              CLAMP((float)(iteration - params.tau)/(params.max_iterations - params.tau), 0, 1);
                    voxels.set_occluder_pyramid2(new_pos, 1, b, params.occlusion_pyramid_d);
                    sp_data.remove_close(new_pos, 0.75*params.ro);
                    
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
        float total_light = 0;
        int total_joints = b.joints.size();
        int i=0;
        for (Joint &j : b.joints)
        {
            if (i > -10)
            {
                total_light += MIN(1, 1/(0.5 + voxels.get_occlusion_simple(j.pos)));
                //logerr("l %f", voxels.get_occlusion_simple(j.pos));
                for (Branch &br : j.childBranches)
                {
                    if (br.alive)
                        remove_branches(t, br, params, voxels);
                    total_light += br.total_resource;//total_light or 0 if branch is removed
                    total_joints += br.total_joints;
                }

/*
                auto it = j.childBranches.begin();
                while (it != j.childBranches.end())
                {
                    if (it->alive)
                    {
                        it++;
                    }
                    else
                    {
                        it = j.childBranches.erase(it);
                    }
                }
*/
            }
            i++;
        }
        float lq = total_light/total_joints;
        if (total_joints > 10)
        {
            //logerr("br is ok %d %d %f",b.level, total_joints, total_light);
        }
        if (b.level > 0 && lq < (params.r_s + 0.03*MIN(log2f(total_joints),10)) && b.joints.size() > 1)
        {
            //remove branch
            //logerr("remove branch %d %f",total_joints, total_light);
            set_occlusion(b, voxels, params, -1);//remove occlusion from deleted branch
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
                r += pow(br.base_r, params.r_pow);//use memory of removed branches too
            }
            it->r = pow(r, 1.0/params.r_pow);
        }
        b.base_r = pow(r, 1.0/params.r_pow);
    }

    void GETreeGenerator::set_occlusion(Branch &b, LightVoxelsCube &voxels, GETreeParameters &params, float mul)
    {
        int i = 0;
        for (Joint &j : b.joints)
        {
            if (i > 0)
            {
                float b = params.b_min + (params.b_max - params.b_min)*
                        CLAMP((float)(j.birth_time - params.tau)/(params.max_iterations - params.tau), 0, 1);

                voxels.set_occluder_pyramid2(j.pos, mul, b, params.occlusion_pyramid_d);
                for (Branch &br : j.childBranches)
                {
                    if (br.alive)
                        set_occlusion(br, voxels, params, mul);
                }
            }
            i++;
        }
    }

    bool GETreeGenerator::SpaceColonizationData::find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos, 
                                                               glm::vec3 dir, float angle,
                                                               glm::vec3 &best_pos, float &best_occ)
    {
        best_occ = 1000;
        float cs = cos(angle);

        std::function<void(glm::vec3 &)> func = [&](glm::vec3 &p){
            if (dot(normalize(p-pos),dir) > cs)
                {
                    float occ = voxels.get_occlusion(p);
                    //logerr("%f occ", occ);
                    if (occ < best_occ)
                    {
                        best_occ = occ;
                        best_pos = p;
                    }
                }
        };

        AABB box = AABB(pos - r*vec3(1,1,1), pos + r*vec3(1,1,1)); 
        octree.apply_to_neighbours_sphere(box, r, pos, func);

        if (best_occ >= 1000)
            return false;
        else 
            return true;
    }

    void GETreeGenerator::SpaceColonizationData::remove_close(glm::vec3 pos, float r)
    {
        AABB box = AABB(pos - r*vec3(1,1,1), pos + r*vec3(1,1,1)); 
        octree.remove_in_sphere(box, r, pos);
    }