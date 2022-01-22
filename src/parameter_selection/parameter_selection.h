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
    }; 
    Results parameter_selection(TreeTypeData reference_tree, Block &selection_settings, Scene *demo_scene = nullptr);
    Results parameter_selection(Block &reference_info, Block &selection_settings, Scene *demo_scene = nullptr);
private:
    void parameter_selection_internal(Block &selection_settings, Results &results, Scene &scene, ReferenceTree &ref_tree,
                                      TreeTypeData *ref_type = nullptr);
};