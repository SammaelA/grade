#include "load_tree_structure.h"
#include <glm/gtc/matrix_transform.hpp>

std::string TreeLoaderBlk::trees_directory = "./original_trees";
std::string TreeLoaderBlk::blks_base_name = "tree";
int TreeLoaderBlk::trees_taken = 0;
int tlb_ids = 0;
int br_ids = 0;
struct SJoint
{
    glm::vec3 pos;
    float r;
    std::vector<int> child_splines;
    std::vector<int> child_leaves;
};
struct SLeaf
{
    glm::vec3 pos;
    glm::vec3 dir;
    glm::vec3 right;
};
struct Spline
{
    std::vector<SJoint> joints;
    int level = 0;
};

struct SavedTree
{
    int trunk_n  = 0;
    glm::vec3 pos = glm::vec3(0,0,0);
    float scale = 1.0;
    std::vector<Spline> splines;
    std::vector<SLeaf> leaves;
    bool on_heightmap = true;

    int joints_count = 0;
    int joints_transformed = 0;
};

void restore_graph(SavedTree &t)
{
    constexpr int sp_hash_sz = 10;
    constexpr float rel_eps = 1e-4;

    glm::vec3 min_pos = glm::vec3(1e9, 1e9, 1e9);
    glm::vec3 max_pos = glm::vec3(-1e9, -1e9, -1e9);

    for (auto &s : t.splines)
    {
        min_pos = min(min_pos, s.joints.front().pos);
        max_pos = max(max_pos, s.joints.front().pos);
        
        min_pos = min(min_pos, s.joints.back().pos);
        max_pos = max(max_pos, s.joints.back().pos);
    }
    typedef std::list<std::pair<int,int>> hash_bucket;
    std::vector<hash_bucket> hash_table = std::vector<hash_bucket>(sp_hash_sz*sp_hash_sz*sp_hash_sz, hash_bucket());
    //glm::vec3 grid_sz = 
    #define P(a,mn,mx) (CLAMP( (int)(((a-mn)/(mx-mn))*sp_hash_sz),0,sp_hash_sz - 1))
    #define HASH(pos, dx, dy, dz) (sp_hash_sz*sp_hash_sz*CLAMP(P((pos).x,min_pos.x,max_pos.x) + dx,0,sp_hash_sz - 1) + \
                                   sp_hash_sz*CLAMP(P((pos).y,min_pos.y,max_pos.y) + dy,0,sp_hash_sz - 1)+ \
                                   CLAMP(P((pos).z,min_pos.z,max_pos.z) + dz,0,sp_hash_sz - 1))
    
    for (int i=0;i<t.splines.size();i++)
    {
        for (int j=0;j<t.splines[i].joints.size();j++)
        {
            int hash_id = HASH(t.splines[i].joints[j].pos,0,0,0);
            hash_table[hash_id].push_back({i,j});
        }
    }

    for (int i=0;i<t.splines.size();i++)
    {
        auto &s = t.splines[i];
        auto &bp = s.joints.front().pos;
        auto &br = s.joints.front().r;
        bool found = false;
        int radius = 0;
        while (!found && radius <= 1)
        {

            float min_dist = pow(100,radius) * rel_eps * glm::dot(max_pos - min_pos, max_pos - min_pos);
            std::pair<int,int> min_pair = {-1,-1};
            for (int dx = -radius; dx <= radius; dx++)
            {
                for (int dy = -radius; dy <= radius; dy++)
                {
                    for (int dz = -radius; dz <= radius; dz++)
                    {
                        int hash_id = HASH(bp, dx, dy, dz);
                        for (auto &p : hash_table[hash_id])
                        {
                            if (p.first == i || p.second == 0)
                                continue;
                            auto &pos = t.splines[p.first].joints[p.second].pos;
                            auto &r = t.splines[p.first].joints[p.second].r;
                            float d = glm::dot(pos - bp, pos - bp);
                            if (d < min_dist && r >= br)
                            {
                                min_dist = d;
                                min_pair = p;
                            }
                        }
                    }
                }
            }

            if (min_pair.first >= 0 && min_pair.second >= 0)
            {
                //bp = t.splines[min_pair.first].joints[min_pair.second].pos;
                t.splines[min_pair.first].joints[min_pair.second].child_splines.push_back(i);
                found = true;
            }
            radius++;
        }
        /*if (!found && i > 0)
        {
            t.splines[0].joints[1].child_splines.push_back(i);
        }*/
    }

    for (int i=0;i<t.leaves.size();i++)
    {
        auto &bp = t.leaves[i].pos;
        int hash_id = HASH(bp, 0, 0, 0);
        float min_dist = 1e10;
        std::pair<int,int> min_pair = {-1,-1};

        for (auto &p : hash_table[hash_id])
        {
            auto &pos = t.splines[p.first].joints[p.second].pos;
            float dist = dot(pos - bp, pos - bp);
            if (dist < min_dist)
            {
                min_dist = dist;
                min_pair = p;
            }
        }
        if (min_pair.first >= 0 && min_pair.second >= 0)
        {
            t.splines[min_pair.first].joints[min_pair.second].child_leaves.push_back(i);
            t.leaves[i].pos = t.splines[min_pair.first].joints[min_pair.second].pos;
        }
    }
}

SavedTree dummy_tree()
{
    SavedTree st;
    st.splines.emplace_back();
    SJoint j1,j2;
    j1.pos = glm::vec3(0,0,0);
    j1.r = 1;
    j2.pos = glm::vec3(0, 10, 0);
    j2.r = 0.5;
    st.splines.back().joints = {j1, j2};
    st.splines.back().level = 0;
}

void load_tree(Block &blk, SavedTree &st)
{
    st.pos = blk.get_vec3("position", st.pos);
    st.on_heightmap = blk.get_bool("on_heightmap", st.on_heightmap);
    st.scale = blk.get_double("scale", st.scale);
    glm::mat4 transform = glm::rotate(glm::mat4(st.scale),-PI/2,glm::vec3(1,0,0));

    Block *splines_bl = blk.get_block("splines");
    if (!splines_bl)
    {
        logerr("given .bsg file does not contain any splines");
        return;
    }
    else
    {
        int splines_cnt  = splines_bl->size();
        for (int i=0;i<splines_cnt;i++)
        {
            std::vector<float> spline_raw;
            bool is_arr = splines_bl->get_arr(i, spline_raw);
            if (is_arr && spline_raw.size() > 0 && (spline_raw.size() % 4 == 0))
            {
                st.splines.emplace_back();
                for (int j=0;j<spline_raw.size();j+=4)
                {
                    glm::vec3 pos = st.pos + glm::vec3(transform * glm::vec4(spline_raw[j], spline_raw[j+1], spline_raw[j+2], 1));
                    float r = spline_raw[j+3] * st.scale;
                    st.splines.back().joints.push_back(
                    SJoint{pos, r, std::vector<int>()});
                    st.joints_count++;
                }
            }
        }
    }

    std::vector<float> leaves_arr;
    blk.get_arr("leaves",leaves_arr);
    if (leaves_arr.size() > 0)
    {
        for (int i=0;i<leaves_arr.size()-8;i+=9)
        {
            st.leaves.emplace_back();
            st.leaves.back().pos = glm::vec3(leaves_arr[i],leaves_arr[i+1],leaves_arr[i+2]);
            st.leaves.back().dir = glm::vec3(leaves_arr[i+3],leaves_arr[i+4],leaves_arr[i+5]);
            st.leaves.back().right = glm::vec3(leaves_arr[i+6],leaves_arr[i+7],leaves_arr[i+8]);

            st.leaves.back().pos = st.pos + glm::vec3(transform*glm::vec4(st.leaves.back().pos,1));
            st.leaves.back().dir = glm::vec3(transform*glm::vec4(st.leaves.back().dir,0));
            st.leaves.back().right = glm::vec3(transform*glm::vec4(st.leaves.back().right,0));
        }
    }
    else
    {
        logerr("tree has malformed leaves array");
    }

    Block *struct_bl = blk.get_block("structure");
    if (struct_bl)
    {
        //TODO
    }
    else
    {
        restore_graph(st);
    }
}

void convert(Spline &s, ::Branch *br, int level, SavedTree &st, ::Tree &external_tree, glm::vec3 &parent_joint_pos)
{
    #define EPS 1e-4
    br->level = level;
    br->dead = false;
    br->self_id = br_ids;
    br_ids++;
    br->id = external_tree.id;
    br->center_self = s.joints[0].pos;
    br->plane_coef = glm::vec4(1,0,0,br->center_self.x);
    if (glm::length(parent_joint_pos - s.joints[0].pos) > 0.5*s.joints[0].r + EPS)
    {
        br->joints.emplace_back();
        br->joints.back().leaf = nullptr;
        br->joints.back().pos = parent_joint_pos;   

        br->segments.emplace_back();
        br->segments.back().begin = parent_joint_pos;
        br->segments.back().rel_r_begin = MAX(EPS,s.joints[0].r);

        br->segments.back().end = s.joints[0].pos;
        br->segments.back().rel_r_end = MAX(EPS,s.joints[0].r);    
    }
    for (int i=0;i<s.joints.size();i++)
    {
        br->joints.emplace_back();
        br->joints.back().leaf = nullptr;
        br->joints.back().pos = s.joints[i].pos;
        st.joints_transformed++;
        //logerr("added joint %f %f %f %f ",s.joints[i].pos.x, s.joints[i].pos.y, s.joints[i].pos.z, s.joints[i].r);
        if (i != s.joints.size()-1)
        {
            br->segments.emplace_back();
            br->segments.back().begin = s.joints[i].pos;
            br->segments.back().rel_r_begin = MAX(EPS,s.joints[i].r);

            br->segments.back().end = s.joints[i+1].pos;
            br->segments.back().rel_r_end = MAX(EPS,s.joints[i+1].r);
        }

        for (int child_spline : s.joints[i].child_splines)
        {
            while (level + 1 >= external_tree.branchHeaps.size())
            {
                external_tree.branchHeaps.push_back(new BranchHeap());
            }
            Branch *ch_b = external_tree.branchHeaps[level + 1]->new_branch();
            ch_b->center_par = br->center_self;
            convert(st.splines[child_spline], ch_b, level+1, st, external_tree, s.joints[i].pos);
            br->joints.back().childBranches.push_back(ch_b);
        }
        
        for (int child_leaf : s.joints[i].child_leaves)
        {
            if (br->joints.back().leaf == nullptr)
            {
                auto leaf = (child_leaf >= 0 && child_leaf < st.leaves.size()) ? st.leaves[child_leaf] : SLeaf();
                br->joints.back().leaf = external_tree.leaves->new_leaf();
                br->joints.back().leaf->pos = leaf.pos + 0.5f*leaf.right;
                br->joints.back().leaf->edges.push_back(leaf.pos);
                br->joints.back().leaf->edges.push_back(leaf.pos + leaf.right);
                br->joints.back().leaf->edges.push_back(leaf.pos + leaf.dir);
                br->joints.back().leaf->edges.push_back(leaf.pos + leaf.dir + leaf.right);
            }
        }
    }
}
void convert_tree(SavedTree &st, GroveGenerationData ggd, ::Tree &external_tree, Heightmap &h)
{
    external_tree.pos = st.pos;
    if (st.on_heightmap)
    {
        external_tree.pos.y = h.get_height(st.pos);
    }
    external_tree.id = tlb_ids;
    external_tree.valid = true;
    external_tree.type = &(ggd.types[0]);
    external_tree.leaves = new LeafHeap();
    external_tree.branchHeaps.push_back(new BranchHeap());
    Branch *root = external_tree.branchHeaps[0]->new_branch();
    root->center_par = external_tree.pos;
    external_tree.root = root;
    convert(st.splines[st.trunk_n], root, 0, st, external_tree, root->center_par);

    tlb_ids++;
}
void TreeLoaderBlk::load_tree(GroveGenerationData ggd, ::Tree &tree_external, Heightmap &h, Block &b)
{
    SavedTree st;
    ::load_tree(b, st);
    convert_tree(st, ggd, tree_external, h);
    logerr("converted tree %d joints loaded %d joint transformed",st.joints_count,st.joints_transformed);
}
void TreeLoaderBlk::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    BlkManager man;
    for (int i=0;i<ggd.trees_count;i++)
    {
        std::string blk_path = trees_directory + "/" + blks_base_name + "_" + std::to_string(trees_taken) + ".bsg";
        Block b;
        man.load_block_from_file(blk_path, b);
        load_tree(ggd,trees_external[i],h,b);
        trees_taken++;
        
        /*for (auto *bh : trees_external[i].branchHeaps)
        {
            for (auto &b : bh->branches)
            {
                for (auto &s : b.segments)
                {
                    logerr("%d %d segment %f %f %f --> %f %f %f",b.level, b.self_id, s.begin.x, s.begin.y, s.begin.z,
                    s.end.x, s.end.y, s.end.z);
                }
            }
        }*/
    }
}