#pragma once
#include "core/tree.h"
#include "core/scene.h"
#include "common_utils/parameter.h"

struct ReferenceTree;
class ParameterSelector
{
public:
    struct Results
    {
        std::vector<TreeTypeData> best_candidates;
        float best_res = 0;
    }; 
    Results parameter_selection(Block &reference_info, Block &selection_settings);
private:
    void parameter_selection_internal(Block &selection_settings, Results &results, ReferenceTree &ref_tree,
                                      TreeTypeData *ref_type = nullptr, bool save_result_image = false);
};