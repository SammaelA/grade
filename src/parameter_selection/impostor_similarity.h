#pragma once
#include "core/scene.h"
#include "graphics_utils/texture_manager.h"
#include "common_utils/utility.h"


enum TCIFeatureStatus
{
    UNKNOWN,
    EXPLICIT,
    DONT_CARE,
    FROM_TYPE,
    FROM_IMAGE
};

struct TreeCompareInfo
{
    glm::vec2 BCyl_sizes = glm::vec2(0,0);//radius and height of bounding cylinder;
    float branches_density = 0;
    float leaves_density = 0;
    float branches_curvature = 0;//from 0 to 1, average dot(seg, seg_next) for all branches
    float trunk_thickness = 0; 
    int joints_cnt = 0;
    float tropism = 0;
};
struct ReferenceTree
{
    TreeCompareInfo info;
    TCIFeatureStatus reference_image_status = EXPLICIT;
    TCIFeatureStatus width_status = EXPLICIT;
    TCIFeatureStatus height_status = EXPLICIT;
    TCIFeatureStatus branches_density_status = EXPLICIT;
    TCIFeatureStatus leaves_density_status = EXPLICIT;
    TCIFeatureStatus branches_curvature_status = EXPLICIT;
    TCIFeatureStatus trunk_thickness_status = EXPLICIT;
    TCIFeatureStatus joints_cnt_status = EXPLICIT;
    TCIFeatureStatus tropism_status = EXPLICIT;
    TextureAtlas atlas;
    TreeTypeData *reference_type = nullptr;
    ReferenceTree(){};
};

class ImpostorSimilarityCalc
{
public:
    ImpostorSimilarityCalc(int max_impostors, int slices_per_impostor, bool use_top_slice = false);
    ~ImpostorSimilarityCalc();
    void calc_similarity(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results,
                         Tree *original_trees, bool debug_print = false, bool image_debug = false);
    static void get_tree_compare_info(Impostor &imp, Tree &original_tree, TreeCompareInfo &info);
private:
    GLuint fbo=0, slices_info_buf=0, results_buf=0, impostors_info_buf=0, dbg_buf=0;
    int slices_per_impostor = 8;
    int slices_stride = 9;
    bool use_top_slice = false;
    int max_impostors = 1024;
    float *results_data = nullptr;
    glm::uvec4 *slices_info_data = nullptr;
    Shader similarity_shader;
    TreeCompareInfo *impostors_info_data = nullptr;
};