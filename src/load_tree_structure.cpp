#include "load_tree_structure.h"

std::string TreeLoaderBlk::trees_directory = "./original_trees";
std::string TreeLoaderBlk::blks_base_name = "tree";
int TreeLoaderBlk::trees_taken = 0;
int tlb_ids = 0;
int br_ids = 0;
struct SJoint
{
    glm::vec3 pos;
    float r;
    int child_spline = -1;
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
    bool on_heightmap = true;
};

void restore_graph(SavedTree &t)
{
    constexpr int sp_hash_sz = 5;
    constexpr float rel_eps = 1e-6;

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
    
    #define P(a,mn,mx) (CLAMP( (int)(((a-mn)/(mx-mn))*sp_hash_sz),0,sp_hash_sz - 1))
    #define HASH(pos) (sp_hash_sz*sp_hash_sz*P((pos).x,min_pos.x,max_pos.x)+sp_hash_sz*P((pos).y,min_pos.y,max_pos.y)+P((pos).z,min_pos.z,max_pos.z))
    
    for (int i=0;i<t.splines.size();i++)
    {
        for (int j=0;j<t.splines[i].joints.size();j++)
        {
            int hash_id = HASH(t.splines[i].joints[j].pos);
            //TODO
            hash_table[hash_id].push_back({i,j});
        }
    }

    for (int i=0;i<t.splines.size();i++)
    {
        auto &s = t.splines[i];
        auto &bp = s.joints.front().pos;
        int hash_id = HASH(bp);

        float min_dist = rel_eps * glm::dot(max_pos - min_pos, max_pos - min_pos);
        std::pair<int,int> min_pair = {-1,-1};
        for (auto &p : hash_table[hash_id])
        {
            if (p.first == i)
                continue;
            auto &pos = t.splines[p.first].joints[p.second].pos;
            float d = glm::dot(pos - bp, pos - bp);
            if (d < min_dist)
            {
                min_dist = d;
                min_pair = p;
            }
        }
        if (min_pair.first >= 0 && min_pair.second >= 0)
        {
            t.splines[min_pair.first].joints[min_pair.second].child_spline = i;
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
                    st.splines.back().joints.push_back(
                    SJoint{glm::vec3(spline_raw[j], spline_raw[j+1], spline_raw[j+2]), spline_raw[j+3], -1});
                }
            }
        }
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

void convert(Spline &s, ::Branch *br, int level, SavedTree &st, ::Tree &external_tree)
{
    br->level = level;
    br->dead = false;
    br->self_id = br_ids;
    br_ids++;
    br->id = external_tree.id;
    br->center_self = s.joints[0].pos;
    br->plane_coef = glm::vec4(1,0,0,br->center_self.x);
    for (int i=0;i<s.joints.size();i++)
    {
        br->joints.emplace_back();
        br->joints.back().leaf = nullptr;
        br->joints.back().pos = s.joints[i].pos;
        //logerr("added joint %f %f %f %f ",s.joints[i].pos.x, s.joints[i].pos.y, s.joints[i].pos.z, s.joints[i].r);
        if (i != s.joints.size()-1)
        {
            br->segments.emplace_back();
            br->segments.back().begin = s.joints[i].pos;
            br->segments.back().rel_r_begin = s.joints[i].r;

            br->segments.back().end = s.joints[i+1].pos;
            br->segments.back().rel_r_end = s.joints[i+1].r;
        }

        if (s.joints[i].child_spline >= 0)
        {
            while (level + 1 >= external_tree.branchHeaps.size())
            {
                external_tree.branchHeaps.push_back(new BranchHeap());
            }
            Branch *ch_b = external_tree.branchHeaps[level + 1]->new_branch();
            ch_b->center_par = br->center_self;
            convert(st.splines[s.joints[i].child_spline], ch_b, level+1, st, external_tree);
            br->joints.back().childBranches.push_back(ch_b);
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
    external_tree.pos = glm::vec3(0,0,0);
    external_tree.id = tlb_ids;
    external_tree.valid = true;
    external_tree.type = &(ggd.types[0]);
    external_tree.leaves = new LeafHeap();
    external_tree.branchHeaps.push_back(new BranchHeap());
    Branch *root = external_tree.branchHeaps[0]->new_branch();
    root->center_par = glm::vec3(0,0,0);
    external_tree.root = root;
    convert(st.splines[st.trunk_n], root, 0, st, external_tree);

    tlb_ids++;
}

void TreeLoaderBlk::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    BlkManager man;
    for (int i=0;i<ggd.trees_count;i++)
    {
        std::string blk_path = trees_directory + "/" + blks_base_name + "_" + std::to_string(trees_taken) + ".bsg";
        Block b;
        SavedTree st;
        man.load_block_from_file(blk_path, b);
        load_tree(b, st);
        convert_tree(st, ggd, trees_external[i], h);
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