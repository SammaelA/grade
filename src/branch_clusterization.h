#pragma once
#include "tree.h"
#include "visualizer.h"
#include <vector>
#include <map>
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
        Answer(): Answer(false,0,1){};
    };
    struct DistDataTable
    {
        Answer *data = nullptr;
        int n;
        void create(int _n)
        {
            n = _n;
            data = new Answer[n*n];
        }
        void clear()
        {
            if (data)
                delete[] data;
        }
        inline Answer get(int x, int y){return data[x*n+y];}
        inline void set(int x, int y, Answer &a){data[x*n+y] = a;}
        int size() {return n;}
        ~DistDataTable() {clear();}
    };
    static void calc_joints_count(Branch *b, std::vector<int>  &counts);
    struct BranchWithData
    {
        Branch *original;
        Branch *b;
        int pos;
        std::vector<int> joint_counts;
        glm::mat4 transform;
        BranchWithData(Branch *_original, Branch *_b = nullptr, int levels = 0, int _pos = 0, glm::mat4 _transform = glm::mat4(1.0f))
        {
            original = _original;
            b = _b;
            pos = _pos;
            transform = _transform;
            for (int i=0;i<levels;i++)
                joint_counts.push_back(0);
            if (b)
                calc_joints_count(b,joint_counts);
        }
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
        BranchWithData *branch = nullptr;
        Cluster *U = nullptr;
        Cluster *V = nullptr;
        int size = 0;
        std::map<Cluster *,DistData> distances;
        Cluster(BranchWithData *bwd);
        Cluster(Cluster *_U, Cluster *_V);
        void to_branch_data(std::vector<Branch *> &branches);
        void to_base_clusters(std::vector<Cluster *> &clusters);
        float ward_dist(Cluster *B, float min = 1.0, float max = 0.0);
        Branch *prepare_to_replace(std::vector<glm::mat4> &transforms);
        
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
            clusters.reserve(branches.size()*2);
            for (int i=0;i<branches.size();i++)
            {
                clusters.push_back(Cluster(&(branches[i])));
                current_clusters.push_back(i);
            }
        }
        void make(int n = 20);
        Dist get_P_delta(int n,std::list<int> &current_clusters, std::list<Dist> &P_delta, float &delta);
    };
    bool set_branches(Tree &t, int layer);
    bool set_branches(Tree *t, int count, int layer);
    void visualize_clusters(DebugVisualizer &debug, bool need_debug = false);
    static bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                             float min, float max);
    static bool match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc, std::vector<int> &jp,
                                     float min, float max);
    static Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0, DistData *data = nullptr);
    static Answer dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    static Answer dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    static Answer dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0, DistData *data = nullptr);
    Answer cluster_dist_min(Cluster &c1, Cluster &c2, float min = 1.0, float max = 0.0);
    Clusterizer()
    {
    }
    static float delta;
    BranchHeap branchHeap;
    LeafHeap leafHeap;
    std::vector<BranchWithData> branches;
    DistDataTable ddt;
    static std::vector<float> weights;
    ClusterDendrogramm Ddg;
}; 