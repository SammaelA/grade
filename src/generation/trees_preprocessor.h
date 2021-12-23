#pragma once
#include "core/tree.h"

class TreePreprocessor
{
public:
    void preprocess_tree(Tree &t, Block &preprocessing_params);
private:
    void simplify_branch_rec(Tree &t, Branch *b, int max_merge_n);
};