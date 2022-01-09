#pragma once
#include "core/scene.h"
#include "graphics_utils/texture_manager.h"
#include "common_utils/utility.h"

struct TreeCompareInfo
{
    glm::vec2 BCyl_sizes;//radius and height of bounding cylinder;
    glm::vec2 _pad;
};
struct ReferenceTree
{
    TreeCompareInfo info;
    Texture tex;
    ReferenceTree() : tex(textureManager.empty()) {};
};

class ImpostorSimilarityCalc
{
public:
    ImpostorSimilarityCalc(int max_impostors, int slices_per_impostor, bool use_top_slice = false);
    ~ImpostorSimilarityCalc();
    void calc_similarity(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results);
    static void get_tree_compare_info(Impostor &imp, TreeCompareInfo &info);
private:
    GLuint fbo, slices_info_buf, results_buf, impostors_info_buf;
    int slices_per_impostor = 8;
    int slices_stride = 9;
    bool use_top_slice = false;
    int max_impostors = 1024;
    float *results_data = nullptr;
    glm::uvec4 *slices_info_data = nullptr;
    Shader similarity_shader;
    TreeCompareInfo *impostors_info_data = nullptr;
};