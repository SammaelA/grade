#include "trees_preprocessor.h"

using namespace glm;

int prev_count = 0;
int new_count = 0;
int t1 = 0;
int t2 = 0;
int t3 = 0;
int t4 = 0;

void TreePreprocessor::preprocess_tree(Tree &t, Block &preprocessing_params)
{
    int max_merge_n = preprocessing_params.get_int("max_merge_n",12);
    simplify_branch_rec(t, t.root, max_merge_n);
    t.leaves->clear_removed();
    for (auto &bh : t.branchHeaps)
        bh->clear_removed();
    logerr("simplified tree %d/%d joints stat %d %d %d %d", new_count, prev_count, t1, t2, t3, t4);
}

void TreePreprocessor::simplify_branch_rec(Tree &t, Branch *b, int max_merge_n)
{
    float r_diff_thr = 0.67;
    float dir_diff_thr = 0.5;
    if (b->dead || b->joints.size() < 2)
        return;
    bool can_be_simplified = true;
    for (auto &s : b->segments)
    {
        if (!s.mults.empty())
        {
            can_be_simplified = false;
            break;
        }
    }
    if (can_be_simplified)
    {
        std::vector<std::list<Joint>::iterator> all_joints;
        std::vector<float> all_rs;
        std::vector<bool> is_useful;

        int merge_n = 0;
        auto it = b->joints.begin();
        auto sit = b->segments.begin();

        all_joints.push_back(it);
        all_rs.push_back(sit->rel_r_begin);
        is_useful.push_back(true);
        auto next_sit = sit;
        next_sit++;
        it++;
        float last_r = sit->rel_r_begin;
        vec3 last_dir = normalize(sit->end - sit->begin);

        while (it != b->joints.end() && sit != b->segments.end())
        {
            bool useful = false;
            vec3 dir = last_dir;
            if (next_sit != b->segments.end())
                dir = normalize(next_sit->end - sit->begin);
            if (!it->childBranches.empty()) //has child branches
            {
                useful = true;
                t1++;
            }
            else if ((next_sit != b->segments.end()) && 
                    ((next_sit->rel_r_end/last_r) < r_diff_thr ||
                      dot(dir, last_dir) < dir_diff_thr))
                    //branch r changes significantly
            {
                useful = true;
                t2++;
            }
            else if (merge_n >= max_merge_n)//to far from previous useful node
            {
                useful = true;
                t3++;
            }
            else if (next_sit == b->segments.end())//end of branch
            {
                useful = true;
                t4++;
            }
            if (useful)
            {
                last_r = next_sit->rel_r_begin;
                last_dir = dir;
            }
            all_joints.push_back(it);
            all_rs.push_back(sit->rel_r_end);
            is_useful.push_back(useful);
            if (useful)
                merge_n = 0;
            else
                merge_n++;
            it++;
            sit++;
            if (next_sit != b->segments.end())
                next_sit++;
        }

        std::list<Joint> new_joints;
        std::list<Segment> new_segments;
        float prev_r = -1;
        for (int i=0;i<is_useful.size();i++)
        {
            if (is_useful[i])
            {
                if (!new_joints.empty())
                {
                    new_segments.emplace_back();
                    new_segments.back().begin = new_joints.back().pos;
                    new_segments.back().end = all_joints[i]->pos;
                    new_segments.back().rel_r_begin = prev_r;
                    new_segments.back().rel_r_end = all_rs[i];
                    new_segments.back().mults = {};
                }

                new_joints.push_back(Joint(*(all_joints[i])));
                prev_r = all_rs[i];
            }
            else
            {
                if (all_joints[i]->leaf)
                {
                    if (!new_joints.back().leaf)
                    {
                        new_joints.back().leaf = all_joints[i]->leaf;
                    }
                    else
                    {
                        for (auto &e : all_joints[i]->leaf->edges)
                            new_joints.back().leaf->edges.push_back(e);
                        all_joints[i]->leaf->edges.clear();
                    }
                }
            }
        }
        prev_count += b->joints.size();
        new_count += new_joints.size();
        b->joints = new_joints;
        b->segments = new_segments;
    }

    for (auto &j : b->joints)
    {
        for (auto *chb : j.childBranches)
            simplify_branch_rec(t, chb, max_merge_n);
    }

}