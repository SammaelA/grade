#pragma once
#include "core/tree.h"

class TreePreprocessor
{
public:
    void preprocess_tree(Tree &t, Block &preprocessing_params);
private:
    void simplify_branch_rec(Tree &t, Branch *b, int max_merge_n);
    int prev_count = 0;
    int new_count = 0;
    int t1 = 0;
    int t2 = 0;
    int t3 = 0;
    int t4 = 0;
};