#include "simpliest_generator.h"
#include "common_utils/distribution.h"

using namespace glm;

std::atomic<int> branch_next_id(0), tree_next_id(0); 

void SimpliestTreeGenerator::plant_tree(glm::vec3 pos, TreeTypeData *type)
{
    tree_positions.push_back(pos);
    types.push_back(type);
}

void SimpliestTreeGenerator::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
    SimpliestTreeStructureParameters dummy_params;
    for (int i=0;i<tree_positions.size();i++)
    {
        vec3 pos = tree_positions[i];
        SimpliestTreeStructureParameters *params = dynamic_cast<SimpliestTreeStructureParameters *>(types[i]->params);
        if (!params)
        {
            logerr("simpliest tree generator got wrong tree type");
            params = &dummy_params;
        }

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

void SimpliestTreeGenerator::create_tree(Tree *tree, vec3 pos, SimpliestTreeStructureParameters &params)
{
    tree->root = tree->branchHeaps[0]->new_branch();
    tree->root->type_id = 0;
    tree->root->self_id = branch_next_id.fetch_add(1);
    tree->root->level = 0;
    tree->root->dead = false;
    tree->root->center_self = pos;
    tree->root->center_par = vec3(0,0,0);
    tree->root->plane_coef = vec4(1,0,0,-pos.x);
    tree->root->id = tree->id;
    create_branch(tree, tree->root, pos, vec3(0,1,0), vec3(1,0,0), 0, params);
    //logerr("created simpliest tree with %d leaves", tree->leaves->leaves.size());
}

void SimpliestTreeGenerator::create_branch(Tree *tree, Branch *branch, glm::vec3 start_pos, glm::vec3 base_dir, 
                                           glm::vec3 normal, int level, 
                                           SimpliestTreeStructureParameters &params)
{
    float seg_len = params.branch_len[level]/params.branch_count[level];
    float r0 = params.branch_r[level];
    float r1 = (level >= params.max_depth - 1) ? r0 : params.branch_r[level + 1];
    vec2 dir_sum = vec2(0,0);
    for (int i=0;i<=params.branch_count[level];i++)
    {
        vec3 pos = start_pos + i*seg_len*base_dir;
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
            if (urand() < params.leaves_count)
            {
                //create leaf
                Leaf *l = tree->leaves->new_leaf();
                l->pos = j.pos;
                glm::vec3 rd1 = rand_dir();
                glm::vec3 rd2 = rand_dir();
                float sz = params.leaf_size;
                glm::vec3 a = j.pos + sz * rd1 + 0.5f * sz * rd2;
                glm::vec3 b = j.pos + 0.5f * sz * rd2;
                glm::vec3 c = j.pos - 0.5f * sz * rd2;
                glm::vec3 d = j.pos + sz * rd1 - 0.5f * sz * rd2;
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
                ch_b->plane_coef = vec4(1,0,0,-j.pos.x);
                ch_b->id = tree->id;
                glm::vec3 nb_dir = base_dir;
                glm::vec3 nb_norm = normal;
                if (i != params.branch_count[level])
                {
                    float psi = params.branch_angle[level];
                    float phi = 2*PI*urand();

                    vec3 z_axis = base_dir;
                    vec3 y_axis = normalize(cross(z_axis, normal));
                    vec3 x_axis = normalize(cross(y_axis, z_axis));

                    float x = sin(psi)*sin(phi);
                    float y = sin(psi)*cos(phi);
                    float z = cos(psi);
                    
                    nb_dir = normalize(x*x_axis + y*y_axis + z*z_axis);
                    float f = length(vec2(nb_dir.x, nb_dir.z));
                    vec2 p_dir = f*normalize(vec2(nb_dir.x, nb_dir.z) - dir_sum);
                    nb_dir = normalize(vec3(p_dir.x, nb_dir.y, p_dir.y));
                    dir_sum += vec2(nb_dir.x, nb_dir.z);
                    //if (level == 0)
                    //logerr("nb %f %f",nb_dir.x, nb_dir.z);
                    nb_norm = normalize(cross(base_dir, nb_dir));
                }
                create_branch(tree, ch_b, j.pos, nb_dir, nb_norm, level+1, params);
                params.set_state(level);
            }
    }
}

void SimpliestTreeStructureParameters::save_to_blk(Block &b)
{
    b.set_arr("branch_len", branch_len);
    b.set_arr("branch_r", branch_r);
    b.set_arr("branch_angle", branch_angle);
    b.set_arr("branch_count", branch_count);
    b.set_arr("branching_start", branching_start);

    b.set_int("max_depth",max_depth);
    b.set_double("leaves_count",leaves_count);
    b.set_double("leaf_size",leaf_size);
}
void SimpliestTreeStructureParameters::load_from_blk(Block &b)
{
    b.get_arr("branch_len", branch_len, true);
    b.get_arr("branch_r", branch_r, true);
    b.get_arr("branch_angle", branch_angle, true);
    b.get_arr("branch_count", branch_count, true);
    b.get_arr("branching_start", branching_start, true);

    max_depth = b.get_int("max_depth",max_depth);
    leaves_count = b.get_double("leaves_count",leaves_count);
    leaf_size = b.get_double("leaf_size",leaf_size);
}

void SimpliestTreeStructureParameters::write_parameter_list(ParameterList &list)
{
    /*    int max_depth = 4;
    std::vector<float> branch_len = {60,20,12,7};
    std::vector<float> branch_r = {1.75,0.75,0.2, 0.075};
    std::vector<float> branch_angle = {PI/6, PI/6, PI/4, PI/4};
    std::vector<int> branch_count = {20,10,8,10};
    std::vector<float> branching_start = {0.5,0.0,0,0};
    float leaves_count = 0.75;
    float leaf_size = 3;*/
    list.ordinalParameters.emplace("max_depth",max_depth);

    list.continuousParameters.emplace("branch_len_0",branch_len[0]);
    list.continuousParameters.emplace("branch_len_1",branch_len[1]);
    list.continuousParameters.emplace("branch_len_2",branch_len[2]);
    list.continuousParameters.emplace("branch_len_3",branch_len[3]);

    list.continuousParameters.emplace("branch_r_0",branch_r[0]);
    list.continuousParameters.emplace("branch_r_1",branch_r[1]);
    list.continuousParameters.emplace("branch_r_2",branch_r[2]);
    list.continuousParameters.emplace("branch_r_3",branch_r[3]);

    list.continuousParameters.emplace("branch_angle_0",branch_angle[0]);
    list.continuousParameters.emplace("branch_angle_1",branch_angle[1]);
    list.continuousParameters.emplace("branch_angle_2",branch_angle[2]);
    list.continuousParameters.emplace("branch_angle_3",branch_angle[3]);

    list.continuousParameters.emplace("branch_count_0",branch_count[0]);
    list.continuousParameters.emplace("branch_count_1",branch_count[1]);
    list.continuousParameters.emplace("branch_count_2",branch_count[2]);
    list.continuousParameters.emplace("branch_count_3",branch_count[3]);

    list.continuousParameters.emplace("branching_start_0",branching_start[0]);
    list.continuousParameters.emplace("branching_start_1",branching_start[1]);
    list.continuousParameters.emplace("branching_start_2",branching_start[2]);
    list.continuousParameters.emplace("branching_start_3",branching_start[3]);

    list.continuousParameters.emplace("leaves_count",leaves_count);
    list.continuousParameters.emplace("leaf_size",leaf_size);

    list.print();
}

void SimpliestTreeStructureParameters::read_parameter_list(ParameterList &list)
{
    max_depth = list.ordinalParameters.at("max_depth");

    branch_len[0] = list.continuousParameters.at("branch_len_0");
    branch_len[1] = list.continuousParameters.at("branch_len_1");
    branch_len[2] = list.continuousParameters.at("branch_len_2");
    branch_len[3] = list.continuousParameters.at("branch_len_3");
    
    branch_r[0] = list.continuousParameters.at("branch_r_0");
    branch_r[1] = list.continuousParameters.at("branch_r_1");
    branch_r[2] = list.continuousParameters.at("branch_r_2");
    branch_r[3] = list.continuousParameters.at("branch_r_3");

    branch_angle[0] = list.continuousParameters.at("branch_angle_0");
    branch_angle[1] = list.continuousParameters.at("branch_angle_1");
    branch_angle[2] = list.continuousParameters.at("branch_angle_2");
    branch_angle[3] = list.continuousParameters.at("branch_angle_3");

    branch_count[0] = list.continuousParameters.at("branch_count_0");
    branch_count[1] = list.continuousParameters.at("branch_count_1");
    branch_count[2] = list.continuousParameters.at("branch_count_2");
    branch_count[3] = list.continuousParameters.at("branch_count_3");

    branching_start[0] = list.continuousParameters.at("branching_start_0");
    branching_start[1] = list.continuousParameters.at("branching_start_1");
    branching_start[2] = list.continuousParameters.at("branching_start_2");
    branching_start[3] = list.continuousParameters.at("branching_start_3");

    leaves_count = list.continuousParameters.at("leaves_count");
    leaf_size = list.continuousParameters.at("leaf_size");
}