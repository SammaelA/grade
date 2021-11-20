#include "trees_preprocessor.h"

void TreePreprocessor::preprocess_tree(Tree &t, Block &preprocessing_params)
{
    int max_merge_n = preprocessing_params.get_int("max_merge_n",5);
    simplify_branch_rec(t.root, max_merge_n);
}

void TreePreprocessor::simplify_branch_rec(Branch *b, int max_merge_n)
{
    float r_diff_thr = 0.75;
    if (b->dead)
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
        while (it != b->joints.end() && sit != b->segments.end())
        {
            bool useful = false;
            if (!it->childBranches.empty() || it->leaf) //has child branches
            {
                useful = true;
            }
            else if ((next_sit != b->segments.end()) && 
                    (next_sit->rel_r_end/next_sit->rel_r_begin) < r_diff_thr)
                    //branch r changes significantly
            {
                useful = true;
            }
            else if (merge_n >= max_merge_n)//to far from previous useful node
            {
                useful = true;
            }
            else if (next_sit == b->segments.end())//end of branch
            {
                useful = true;
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
        }
        b->joints = new_joints;
        b->segments = new_segments;
    }

    for (auto &j : b->joints)
    {
        for (auto *chb : j.childBranches)
            simplify_branch_rec(chb,max_merge_n);
    }

}