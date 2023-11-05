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
    static void visualize_tree(const TreeTypeData &tree_type, const std::string &file_name,
                               int image_count, float distance, glm::ivec2 image_size, int rays_per_pixel,
                               bool render_terrain);
    static void create_single_tree_scene(const TreeTypeData &tree_type, Scene &outScene);
    static void save_tree_to_obj(const TreeTypeData &tree_type, const std::string &file_name);
private:
    void parameter_selection_internal(Block &selection_settings, Results &results, ReferenceTree &ref_tree,
                                      TreeTypeData *ref_type = nullptr, bool save_result_image = false);
};