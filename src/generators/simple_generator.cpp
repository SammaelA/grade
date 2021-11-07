#include "simple_generator.h"

using namespace glm;

vec3 rand_dir()
{
    vec3 dir;
    dir.x = urand(-1,1);
    dir.y = urand(-1,1);
    dir.z = urand(-1,1);
    return normalize(dir);
}

int ids = 0;
int t_ids = 0;

void SimpleTreeGenerator::create_branch(Tree *tree, Branch *branch, vec3 start_pos, vec3 base_dir, vec3 normal, int level, 
                                        float base_r, float leaves_chance)
{

    params.set_state(level);

    int seg_count = params().segment_count();
    //more diversity to test clustering 
    if (level == 1)
    {
        float lrnd = urand();
        leaves_chance = lrnd > 0.4 ? lrnd : lrnd;
        float rnd = urand();
        seg_count = rnd < 0.4 ? (rnd*seg_count + 2) : seg_count;
    }
    branch->joints.emplace_back();
    branch->joints.back().pos = start_pos;
    base_r = 0.67*MIN(1, seg_count/params().segment_count.get_base())*base_r;
    float prev_r = base_r;
    int prev_sector = -1;
    for (int i=0;i<seg_count;i++)
    {
        float len = params().segment_size();
        len *= (i == 0 && level < 2) ? 3 : 1;
        float k = (float)(i+1)/(MAX(params().segment_count.get_base(),seg_count) + 1);
        float k1 = (float)(i+1)/seg_count;
        float r = base_r*(1-k1);
        r = MAX(r, 0.01);
        float a1 = params().base_dir_mult();
        float a2 = params().rand_dir_mult();
        float a3 = params().up_dir_mult()*pow(1 - k,params().up_dir_pow());
        float a4 = params().down_dir_mult()*pow(k, params().down_dir_pow());
        vec3 dir = a1*base_dir + a2*rand_dir() + a3*vec3(0, 1, 0) + a4*vec3(0,-1,0);
        dir = normalize(dir);

        vec3 new_pos = branch->joints.back().pos + len*dir;
        
        branch->segments.emplace_back();
        branch->segments.back().begin = branch->joints.back().pos;
        branch->segments.back().end = new_pos;
        branch->segments.back().rel_r_begin = prev_r;
        branch->segments.back().rel_r_end = r;
        
        prev_r = r;

        branch->joints.emplace_back();
        branch->joints.back().pos = new_pos;

        if (level + 1 == params().max_depth() && urand() < leaves_chance)
        {
            //create leaf
            Joint &j = branch->joints.back();
            Leaf *l = tree->leaves->new_leaf();
            l->pos = j.pos;
            glm::vec3 rd1 = rand_dir();
            glm::vec3 rd2 = rand_dir();
            float sz = params().leaf_mult();
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
        else if (level + 1 < params().max_depth() && (urand() < params().branching_chance()))
        {
            //create branch
            Branch *ch_b = tree->branchHeaps[level + 1]->new_branch();
            branch->joints.back().childBranches.push_back(ch_b);
            ch_b->type_id = branch->type_id;
            ch_b->self_id = ids++;
            ch_b->level = level+1;
            ch_b->dead = false;
            ch_b->center_self = new_pos;
            ch_b->center_par = branch->center_self;
            ch_b->plane_coef = vec4(1,0,0,-new_pos.x);
            ch_b->id = tree->id;

            float psi = params().base_branch_angle();
            int sc = params().sectors_count();
            int sector = round(urand(0, sc - 1));
            if (sector >= prev_sector)
                sector++;
            prev_sector = sector;
            float phi = (2*PI*sector)/sc + urand(0, (2*PI)/sc);

            vec3 z_axis = base_dir;
            vec3 y_axis = normalize(cross(z_axis, normal));
            vec3 x_axis = normalize(cross(y_axis, z_axis));

            float x = sin(psi)*sin(phi);
            float y = sin(psi)*cos(phi);
            float z = cos(psi);

            glm::vec3 nb_dir = normalize(x*x_axis + y*y_axis + z*z_axis);
            glm::vec3 nb_norm = normalize(cross(base_dir, nb_dir));
            create_branch(tree, ch_b, new_pos, nb_dir, nb_norm, level+1, r, leaves_chance);
            params.set_state(level);
        }
    }
}

void SimpleTreeGenerator::create_tree(Tree *tree, vec3 pos)
{
    tree->root = tree->branchHeaps[0]->new_branch();
    tree->root->type_id = 0;
    tree->root->self_id = ids++;
    tree->root->level = 0;
    tree->root->dead = false;
    tree->root->center_self = pos;
    tree->root->center_par = vec3(0,0,0);
    tree->root->plane_coef = vec4(1,0,0,-pos.x);
    tree->root->id = tree->id;
    create_branch(tree, tree->root, pos, vec3(0,1,0), vec3(1,0,0), 0, params().base_thickness(), 0);
}

void SimpleTreeGenerator::create_grove(GroveGenerationData _ggd, ::Tree *trees_external, Heightmap &_h)
{
    ggd = _ggd;
    h = &_h;
    SimpleTreeStructureParameters s_params;
    params = BaseParameterSetWrapper<SimpleTreeStructureParameters>(s_params, s_params.max_depth());
    params.set_state(0);
    for (int i=0;i<ggd.trees_count;i++)
    {
        vec3 pos = vec3(50*(t_ids % 10), 0, 50 * (t_ids / 10));
        pos.y = h->get_height(pos);

        for (int j=0;j<params().max_depth();j++)
        {
            BranchHeap *br = new BranchHeap();
            trees_external[i].branchHeaps.push_back(br);
        }

        trees_external[i].leaves = new LeafHeap();
        trees_external[i].id = t_ids;
        trees_external[i].pos = pos;
        trees_external[i].type = &(ggd.types[0]);
        trees_external[i].valid = true;

        create_tree(trees_external + i, pos);
        t_ids++;
    }
}

void SimpleTreeGenerator::plant_tree(glm::vec3 pos, TreeTypeData *type)
{
    tree_positions.push_back(pos);
    types.push_back(type);
}

void SimpleTreeGenerator::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
    SimpleTreeStructureParameters s_params;
    params = BaseParameterSetWrapper<SimpleTreeStructureParameters>(s_params, s_params.max_depth());
    params.set_state(0);
    for (int i=0;i<tree_positions.size();i++)
    {
        vec3 pos = tree_positions[i];

        for (int j=0;j<params().max_depth();j++)
        {
            BranchHeap *br = new BranchHeap();
            trees_external[i].branchHeaps.push_back(br);
        }

        trees_external[i].leaves = new LeafHeap();
        trees_external[i].id = t_ids;
        trees_external[i].pos = pos;
        trees_external[i].type = types[i];
        trees_external[i].valid = true;

        create_tree(trees_external + i, pos);
        t_ids++;
    }
}