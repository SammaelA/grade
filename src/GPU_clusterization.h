#include "branch_clusterization.h"
#include "limits.h"
#define uint unsigned //CPU
#define uvec4 glm::uvec4 //CPU
#define CNT_bytes 8
#define MAX_JOINTS 42
#define _cnt(x) ((x) & (1<<CNT_bytes))
#define _id(x) ((x) >> CNT_bytes)
#define id_cnt_pack(id, cnt) ((id>>CNT_bytes) | cnt)

class GPUClusterizationHelper
{
public:
    void prepare_ddt(std::vector<Clusterizer::BranchWithData> &branches, Clusterizer::DistDataTable &ddt, ClusterizationParams &cp);
    GPUClusterizationHelper();
    ~GPUClusterizationHelper();

private:
    struct gJoint
    {
        uint r_id;
        uint pos_id;
        uint child_branches_ids[MAX_CHILD_BRANCHES];
    };// = uint4
    struct gBranch
    {
        uint joint_offset;
        uint joint_count;
        float self_weight;
        float cumulative_weight;
    };
    struct gBranchWithData
    {
        uvec4 voxels_xyz_rots;
        
        uint branch_pos;//position of root gBranch in buffer
        uint branch_id;
        uint voxels_offset;
        uint voxels_size;
        
        uint level;
        uint type_id;
        uint dead;
        uint pad1;
        
        uint joints_counts[MAX_BRANCH_LEVELS];
        uint pad[4*(MAX_BRANCH_LEVELS/4 + 1) - MAX_BRANCH_LEVELS];
    };
    struct gClusterizationParams
    {
        int bwd_rotations = 18;
        int r_samples = 4;
        float delta = 0.2;
        float light_importance = 0.4;
        float voxels_size_mult = 1/2.5;
        int ignore_structure_level = 1000;
        int min_clusters = 1;
        float max_individual_dist = 0.95;
        bool different_types_tolerance = true;
        float weights[MAX_BRANCH_LEVELS];
        float light_weights[MAX_BRANCH_LEVELS];
        float r_weights[MAX_BRANCH_LEVELS];
    };
    std::vector<glm::vec4> positions;
    float *all_voxels = nullptr;
    std::vector<float> joint_rs;
    std::vector<gBranch> sticks;
    std::vector<gBranchWithData> branches;
    std::vector<gJoint> joints;
    float *dist_data = nullptr;
    int branches_size;
    int cur_voxels_pointer;
    gClusterizationParams params;
    GLuint pos_buf,voxels_buf,rs_buf,sticks_buf,branches_buf,dist_data_buf,params_buf,
           joints_buf;
    Shader distCompute;
    #define dd_r(i,j) (dist_data[2*(i*branches_size + j) + 1])
    #define dd_dist(i,j) (dist_data[2*(i*branches_size + j)])

    //CPU-only data
    std::vector<glm::mat4> rotates_transforms;

    void fill_branch_data(Clusterizer::BranchWithData &branch);
    uint fill_branch_data(Branch *branch);
    float calc_cumulative_weight(uint stick_id, int level);
    void calculate_distances(int hardness);
    void calculate_dist(int i, int j);
    void setup_buffers();
    glm::vec2 calculate_dist_simple(int i, int j, int rot, float max_dist);
    float calculate_dist_structure(int i, int j, int rot, float max_dist);
    float calculate_dist_light(int i, int j, int rot, float max_dist);
};