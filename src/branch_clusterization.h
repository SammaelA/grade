#pragma once
#include "tree.h"
#include "visualizer.h"
#include "volumetric_occlusion.h"
#include <vector>
#include <map>
#include "tinyEngine/save_utils/blk.h"
#include "hash.h"

struct ClusterizationParams
{
    int bwd_rotations = 18;
    float delta = 0.2;
    float light_importance = 0.4;
    float voxels_size_mult = 1/2.5;
    bool voxelized_structure = false;
    bool hash_dist = false;
    int EV_hasing_cells = 5;
    int EV_hasing_voxels_per_cell = 25; 
    float structure_voxels_size_mult = 1/2.5;
    int ignore_structure_level = 1000;
    int min_clusters = 1;
    float average_cluster_size_goal = 0;
    float max_individual_dist = 0.95;
    bool different_types_tolerance = true;
    std::vector<float> weights = std::vector<float>{5000,800,40,1,0.01};
    std::vector<float> light_weights = std::vector<float>{5000,800,40,1,0.01};
    std::vector<float> r_weights = std::vector<float>{0.5,0.2,0,0,0};
    void load_from_block(Block *b);
};
extern ClusterizationParams clusterizationParams;
struct AdditionalClusterDataArrays
{
    std::vector<float> rotations;
    std::vector<Branch *> originals;
    std::vector<int> ids;
};
struct ClusterData
{
    long id = -1;
    Branch *base = nullptr;
    InstanceDataArrays IDA;
    AdditionalClusterDataArrays ACDA;
    ClusterData();
};
class Clusterizer
{
public:
    struct Answer
    {
        bool exact;
        float from;
        float to;
        Answer(bool ex, float fr, float t)
        {
            exact = ex;
            from = fr;
            to = t;
        }
        Answer() : Answer(false, 0, 1){};
        Answer(const Answer&) = default;
        Answer(Answer&&) = default;
        Answer& operator=(const Answer&) = default;
        Answer& operator=(Answer&&) = default;
        Answer operator*(const float mult)
        {
            return Answer(exact, MIN(from*mult, to*mult), MAX(from*mult, to*mult));
        }
        Answer operator+(const Answer &add)
        {
            return Answer(exact && add.exact, from + add.from, to + add.to);
        }
        Answer operator-(const Answer &sub)
        {
            return Answer(exact && sub.exact, from - sub.to, to -sub.from);
        }
    };

    static void calc_joints_count(Branch *b, std::vector<int> &counts);
    struct BranchWithData
    {
        Branch *original;
        Branch *b;
        int base_cluster_id;
        int id;
        float rot_angle = 0.0;
        std::vector<int> joint_counts;
        glm::mat4 transform;
        std::vector<LightVoxelsCube *> leavesDensity;
        std::vector<LightVoxelsCube *> voxelizedStructures;
        std::vector<Hash> hashes;
        void set_eigen_values_hash(LightVoxelsCube *voxels, Hash &hash, int cells, int voxels_per_cell, int sz);
        void set_occlusion(Branch *b, LightVoxelsCube *light);
        BranchWithData(Branch *_original, Branch *_b, int _base_cluster_id, int levels, int _id, glm::mat4 _transform);
        ~BranchWithData();
        void clear();
    };
    struct DistData
    {
        float dist = 1;
        float rotation = 0;
        DistData(float _dist = 1, float _rotation = 0)
        {
            dist = _dist;
            rotation = _rotation;
        }
    };
    struct Cluster
    {
        static Clusterizer *currentClusterizer; 
        BranchWithData *branch = nullptr;
        Cluster *U = nullptr;
        Cluster *V = nullptr;
        int size = 0;
        std::map<Cluster *, DistData> distances;
        Cluster(BranchWithData *bwd);
        Cluster(Cluster *_U, Cluster *_V);
        void to_branch_data(std::vector<Branch *> &branches);
        void to_base_clusters(std::vector<Cluster *> &clusters);
        float ward_dist(Cluster *B, float min = 1.0, float max = 0.0);
        Branch *prepare_to_replace(std::vector<ClusterData> &base_clusters, InstanceDataArrays &IDA, 
                                   AdditionalClusterDataArrays &ADCA, long &cluster_id);
        Branch *prepare_to_replace(std::vector<ClusterData> &base_clusters, InstanceDataArrays &IDA, 
                                   AdditionalClusterDataArrays &ADCA, std::vector<Cluster *> &clusters, long &cluster_id);
        BranchWithData *get_typical(std::vector<Cluster *> &clusters);
    };
    struct ClusterDendrogramm
    {
        struct Dist
        {
            int U;
            int V;
            float d;
            Dist(int _U, int _V, float _d)
            {
                U = _U;
                V = _V;
                d = _d;
            }
        };
        int size;
        Cluster *root;
        std::vector<Cluster> clusters;
        std::list<int> current_clusters;
        void make_base_clusters(std::vector<BranchWithData> &branches)
        {
            size = branches.size();
            clusters.reserve(branches.size() * 2);
            for (int i = 0; i < branches.size(); i++)
            {
                clusters.push_back(Cluster(&(branches[i])));
                current_clusters.push_back(i);
            }
        }
        void make(int n = 20, int clusters_num = 1);
        Dist get_P_delta(int n, std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta);
    };
    struct DistDataTable
    {
        std::pair<Answer, DistData> *data = nullptr;
        int n;
        void create(int _n)
        {
            n = _n;
            data = safe_new<std::pair<Answer, DistData>>(n * n, "DistDataTable_data");
        }
        void clear()
        {
            safe_delete<std::pair<Answer, DistData>>(data,"DistDataTable_data");
        }
        inline std::pair<Answer, DistData> get(int x, int y) { return data[x * n + y]; }
        inline void set(int x, int y, Answer &a, DistData &d) 
        { data[x * n + y] = std::pair<Answer, DistData>(a,d); }
        int size() { return n; }
        void copy(DistDataTable &ddt)
        {
            n = ddt.n;
            data = safe_new<std::pair<Answer, DistData>>(n * n, "DistDataTable_data");
            memcpy(data, ddt.data, n*n*sizeof(std::pair<Answer, DistData>));
        }
        ~DistDataTable() { clear(); }
    };

    struct ClusterizationTmpData
    {
        std::vector<float> light_weights;
        std::vector<float> weights;
        BranchHeap branchHeap;
        LeafHeap leafHeap;
        std::vector<BranchWithData> branches;
        DistDataTable ddt;
        std::map<int,int> pos_in_table_by_id;
        ClusterDendrogramm Ddg;
    };

    enum ClusteringStep
    {
        TRUNKS,
        BRANCHES,
        TREES,
        CLUSTERING_STEPS_COUNT
    };
    void get_current_ddt(DistDataTable &ddt);
    std::map<int,int> get_id_pos_map();
    static void set_default_clustering_params(ClusterizationParams &params, ClusteringStep step);
    void set_branches(std::vector<ClusterData> &base_clusters);
    void set_light(LightVoxelsCube *_light);
    void get_base_clusters(Tree &t, int layer, std::vector<ClusterData> &base_clusters);
    void get_base_clusters(Tree *t, int count, int layer, std::vector<ClusterData> &base_clusters);
    void visualize_clusters(DebugVisualizer &debug, bool need_debug = false);
    void prepare_ddt();
    void clusterize(ClusterizationParams &params, std::vector<ClusterData> &base_clusters, std::vector<ClusterData> &clusters,
                    bool need_only_ddt = false);
    ClusterData extract_data(std::vector<ClusterData> &base_clusters, Cluster &cl);
    void get_light(Branch *b, std::vector<float> &light, glm::mat4 &transform);
    Answer light_difference(BranchWithData &bwd1, BranchWithData &bwd2);
    bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                             float min, float max);
    bool match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                                     float min, float max);
    Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0, DistData *data = nullptr);
    Answer dist_trunc(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0, DistData *data = nullptr);
    Answer dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    Answer dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    Answer dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0, DistData *data = nullptr);
    Answer cluster_dist_min(Cluster &c1, Cluster &c2, float min = 1.0, float max = 0.0);
    Answer get_dist(BranchWithData &bwd1, BranchWithData &bwd2, DistData *data = nullptr);
    Clusterizer() {};
    ~Clusterizer();
    LightVoxelsCube *current_light = nullptr;
    ClusterizationTmpData *current_data = nullptr;
    static glm::vec3 canonical_bbox;
};
extern LightVoxelsCube *test_voxels_cube[10];