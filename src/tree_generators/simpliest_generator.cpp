#include "simpliest_generator.h"
#include "common_utils/distribution.h"
#include "common_utils/parameter_save_load_defines.h"



#define NO_RANDOM_GEN 1

void SimpliestTreeGenerator::plant_tree(float3 pos, const TreeTypeData *type)
{
    tree_positions.push_back(pos);
    types.push_back(type);
}

void SimpliestTreeGenerator::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
    SimpliestTreeStructureParameters dummy_params;
    for (int i=0;i<tree_positions.size();i++)
    {
        float3 pos = tree_positions[i];
        SimpliestTreeStructureParameters *params = dynamic_cast<SimpliestTreeStructureParameters *>(types[i]->get_params());
        if (!params)
        {
            logerr("simpliest tree generator got wrong tree type");
            params = &dummy_params;
        }
        params->max_depth = MIN(params->max_depth, 4);

        for (int j=0;j<params->max_depth;j++)
        {
            BranchHeap *br = new BranchHeap();
            trees_external[i].branchHeaps.push_back(br);
        }

        trees_external[i].leaves = new LeafHeap();
        trees_external[i].id = tree_next_id.fetch_add(1);
        trees_external[i].pos = pos;
        trees_external[i].type = types[i];
        trees_external[i].valid = true;

        create_tree(trees_external + i, pos, *params);
    }
    tree_positions.clear();
    types.clear();
}

void SimpliestTreeGenerator::create_tree(Tree *tree, float3 pos, SimpliestTreeStructureParameters &params)
{
    int leaves_tries = 0;
    int leaves_cnt = 0;
    tree->root = tree->branchHeaps[0]->new_branch();
    tree->root->type_id = tree->type->type_id;
    tree->root->self_id = branch_next_id.fetch_add(1);
    tree->root->level = 0;
    tree->root->dead = false;
    tree->root->center_self = pos;
    tree->root->center_par = float3(0,0,0);
    tree->root->plane_coef = float4(1,0,0,-pos.x);
    tree->root->id = tree->id;
    create_branch(tree, tree->root, pos, float3(0,1,0), float3(1,0,0), 0, params, leaves_tries, leaves_cnt);
    //logerr("created simpliest tree with %d  %f", tree->leaves->leaves.size(), params.leaves_count);
}

void SimpliestTreeGenerator::create_branch(Tree *tree, Branch *branch, float3 start_pos, float3 base_dir, 
                                           float3 normal, int level, 
                                           SimpliestTreeStructureParameters &params, int &leaves_tries, int &leaves_cnt)
{
    float seg_len = params.branch_len[level]/params.branch_count[level];
    float r0 = params.branch_r[level];
    float r1 = (level >= params.max_depth - 1) ? r0 : params.branch_r[level + 1];
    float2 dir_sum = float2(0,0);
    for (int i=0;i<=params.branch_count[level];i++)
    {
        float3 pos = start_pos + i*seg_len*base_dir;
        if (!branch->joints.empty())
        {
            branch->segments.emplace_back();
            branch->segments.back().begin = branch->joints.back().pos;
            branch->segments.back().end = pos;

            branch->segments.back().rel_r_begin = r0*(1 - (float)(i-1)/params.branch_count[level]) + r1*(float)(i-1)/params.branch_count[level];
            branch->segments.back().rel_r_end = r0*(1 - (float)i/params.branch_count[level]) + r1*(float)i/params.branch_count[level];
 
        }
        branch->joints.emplace_back();
        branch->joints.back().pos = pos;
        Joint &j = branch->joints.back();
        if (i == 0)
            continue;
        if (level == params.max_depth - 1)
        {
            #if NO_RANDOM_GEN > 0
                leaves_tries++;
                bool leaf_b = leaves_cnt/(float)leaves_tries < params.leaves_count;
            #else
                bool leaf_b = urand() < params.leaves_count;// random version
            #endif
            if (leaf_b)
            {
                leaves_cnt++;
                //create leaf
                Leaf *l = tree->leaves->new_leaf();
                l->pos = j.pos;
                #if NO_RANDOM_GEN > 0
                    float psi = PI / 2;
                    float phi = (2 * PI * (i % 6)) / 6;
                    float3 z_axis = base_dir;
                    float3 y_axis = normalize(cross(z_axis, normal));
                    float3 x_axis = normalize(cross(y_axis, z_axis));

                    float x = sin(psi) * sin(phi);
                    float y = sin(psi) * cos(phi);
                    float z = cos(psi);

                    float3 rd1 = normalize(sin(psi) * sin(phi) * x_axis + sin(psi) * cos(phi) * y_axis + cos(psi) * z_axis);
                    float3 rd2 = normalize(sin(psi) * cos(phi) * x_axis + sin(psi) * sin(phi) * y_axis + cos(psi) * z_axis);
                #else
                    float3 rd1 = rand_dir();
                    float3 rd2 = rand_dir();
                #endif
                float sz = params.leaf_size;
                float3 a = j.pos + sz * rd1 + 0.5f * sz * rd2;
                float3 b = j.pos + 0.5f * sz * rd2;
                float3 c = j.pos - 0.5f * sz * rd2;
                float3 d = j.pos + sz * rd1 - 0.5f * sz * rd2;
                l->edges.push_back(a);
                l->edges.push_back(b);
                l->edges.push_back(c);
                l->edges.push_back(d);

                j.leaf = l;
            }
        }
            else if ((float)i/params.branch_count[level] >= params.branching_start[level])
            {
                //create branch
                Branch *ch_b = tree->branchHeaps[level + 1]->new_branch();
                branch->joints.back().childBranches.push_back(ch_b);
                ch_b->type_id = branch->type_id;
                ch_b->self_id = branch_next_id.fetch_add(1);
                ch_b->level = level+1;
                ch_b->dead = false;
                ch_b->center_self = j.pos;
                ch_b->center_par = branch->center_self;
                ch_b->plane_coef = float4(1,0,0,-j.pos.x);
                ch_b->id = tree->id;
                float3 nb_dir = base_dir;
                float3 nb_norm = normal;
                if (i != params.branch_count[level])
                {
                    float psi = params.branch_angle[level];
                    
                    #if NO_RANDOM_GEN > 0
                        float phi = (2*PI*(i % 6))/6;//determined version
                    #else
                        float phi = 2*PI*urand();//random version
                    #endif

                    float3 z_axis = base_dir;
                    float3 y_axis = normalize(cross(z_axis, normal));
                    float3 x_axis = normalize(cross(y_axis, z_axis));

                    float x = sin(psi)*sin(phi);
                    float y = sin(psi)*cos(phi);
                    float z = cos(psi);
                    
                    nb_dir = normalize(x*x_axis + y*y_axis + z*z_axis);
                    float f = length(float2(nb_dir.x, nb_dir.z));
                    float2 p_dir = f*normalize(float2(nb_dir.x, nb_dir.z) - dir_sum);
                    nb_dir = normalize(float3(p_dir.x, nb_dir.y, p_dir.y));
                    dir_sum += float2(nb_dir.x, nb_dir.z);
                    nb_norm = normalize(cross(base_dir, nb_dir));
                }
                create_branch(tree, ch_b, j.pos, nb_dir, nb_norm, level+1, params, leaves_tries, leaves_cnt);
            }
    }
}

void SimpliestTreeStructureParameters::save_load_define(SaveLoadMode mode, Block &b, ParameterList &list)
{
  P_INT(max_depth, mode);
  P_VEC4(branch_len, mode);
  P_VEC4(branch_r, mode);
  P_VEC4(branch_angle, mode);
  P_VEC4(branch_count, mode);
  P_VEC4(branching_start, mode);
  P_FLOAT(leaves_count, mode);
  P_FLOAT(leaf_size, mode);
}