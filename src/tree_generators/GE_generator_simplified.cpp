#include "GE_generator_simplified.h"
#include <algorithm>
using namespace glm;

bool GETreeGeneratorSimplified::iterate(LightVoxelsCube &voxels)
{
    bool growing = false;
    for (auto &t : trees)
    {
        GETreeParameters *p_ptr = dynamic_cast<GETreeParameters *>(t.type->get_params());
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
                t.status = TreeStatus::GROWN;
            else
            {
                std::vector<GrowPoint> growth_points;
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
                grow_nodes(t, params, growth_points, voxels, max_growth);
                remove_branches(t, t.root, params, voxels);

                t.iteration++;
                growing = true;
            }
        }
    }

    return growing;
}

void GETreeGeneratorSimplified::prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params,
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
        }
        else if (j.can_have_child_branches && j.childBranches.size() < params.max_branches)
        {
            growth_points.push_back(GrowPoint(&j, &b, GrowthType::BRANCHING, pd, resource,i));
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

bool GETreeGeneratorSimplified::find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos,
                                              glm::vec3 dir, float angle,
                                              glm::vec3 &best_pos, float &best_occ)
{
    best_occ = 1000;
    float cs = cos(angle);

    std::function<void(glm::vec3 &)> func = [&](glm::vec3 &p)
    {
        if (dot(normalize(p - pos), dir) > cs)
        {
            float occ = voxels.get_occlusion_simple(p);
            //logerr("%f occ", occ);
            if (occ < best_occ)
            {
                best_occ = occ;
                best_pos = p;
            }
        }
    };

    AABB box = AABB(pos - r * vec3(1, 1, 1), pos + r * vec3(1, 1, 1));
    for (int i=0;i<4;i++)
    {
        glm::vec3 p = pos + r * glm::vec3(self_rand(-1,1), self_rand(-1,1), self_rand(-1,1));
        func(p);
    }

    if (best_occ >= 1000)
        return false;
    else
        return true;
}