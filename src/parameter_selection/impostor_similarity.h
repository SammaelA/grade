#pragma once
#include "core/scene.h"
#include "tinyEngine/engine.h"
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
    int id = -1;
    glm::vec2 BCyl_sizes = glm::vec2(0,0);//radius and height of bounding cylinder;
    glm::vec2 real_sizes = glm::vec2(0,0);
    float branches_density = 0;
    float leaves_density = 0;
    float branches_curvature = 0;//from 0 to 1, average dot(seg, seg_next) for all branches
    float trunk_thickness = 0; 
    int joints_cnt = 0;
    float tropism = 0;
};

struct TrunkInfo
{
    float average_thickness = 0;
    float symmetry_x = 0;//in texture coordinates i.e. in [0,1]
    float trunk_split_level = 0;//as a share of tree size
    float trunk_len = 0;//as a share of tree size
    float thickness[16];//
};

struct TreeImageInfo
{
    glm::vec4 tc_transform = glm::vec4(0,0,1,1);//shift and scale, tranform impostor to occupy the whole texture
    float crown_start_level = 0;//as a share of tree size
    float trunk_thickness = 0;//as a share of tree size
    float crown_leaves_share = 0;
    float crown_branches_share = 0;
    TrunkInfo trunk_info;
};

struct StripeInfo
{
    glm::vec4 crown_bord;//(x_start, dense crown x_start, dense crown x_end, x_end)
    glm::vec2 crown_ymin_ymax;
    float branches_share;
    float leaves_share;
};

struct ReferenceTree
{
    static constexpr float border_size = 0.0f;
    TreeCompareInfo info;
    TreeImageInfo image_info;
    std::vector<StripeInfo> alternative_image_info;
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
                         Tree *original_trees, int original_trees_cnt, bool debug_print = false, bool image_debug = false);
    static void get_tree_compare_info(Impostor &imp, Tree &original_tree, TreeCompareInfo &info);
    void get_tree_image_info(TextureAtlas &images_atl, std::map<int, TreeImageInfo> &results, bool image_debug = false);
    void get_reference_tree_image_info(ReferenceTree &reference, float clsh_mult = 0.85);
    
    std::vector<std::vector<StripeInfo>> get_alternative_tree_image_info(TextureAtlas &images_atl, const std::vector<std::pair<int,int>> &slice_id_impostor_n, 
                                                                         TreeCompareInfo *impostors_info, ReferenceTree *reference = nullptr,
                                                                         glm::vec4 *tc_transform = nullptr);
    void get_reference_tree_image_info_alt(ReferenceTree &reference, float clsh_mult = 0.85);
    void calc_similarity_alt(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results,
                             Tree *original_trees, int original_trees_cnt, bool debug_print = false, bool image_debug = false);
private:
    void set_slices_data(GrovePacked &grove);
    void ref_atlas_transform(TextureAtlas &atl);
    GLuint fbo=0, slices_info_buf=0, results_buf=0, impostors_info_buf=0, dbg_buf=0, tree_image_info_buf=0, stripes_results_buf = 0,
    stripes_info_buf = 0;
    int slices_per_impostor = 8;
    int slices_stride = 9;
    bool use_top_slice = false;
    int max_impostors = 1024;
    float *results_data = nullptr;
    glm::uvec4 *slices_info_data = nullptr;
    Shader similarity_shader, tree_info_shader;
    TreeCompareInfo *impostors_info_data = nullptr;
    TreeImageInfo *tree_image_info_data = nullptr;
    StripeInfo *stripes_data = nullptr;
    glm::vec4 *shader_imp_data = nullptr;
};