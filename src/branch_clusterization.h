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
        Branch *b;
        int pos;
        std::vector<int> joint_counts;
        BranchWithData(Branch *_b = nullptr, int levels = 0, int pos = 0)
        {
            b = _b;
            for (int i=0;i<levels;i++)
                joint_counts.push_back(0);
            if (b)
                calc_joints_count(b,joint_counts);
        }
    };
    struct Cluster
    {
        BranchWithData *branch = nullptr;
        Cluster *U = nullptr;
        Cluster *V = nullptr;
        int size = 0;
        std::map<Cluster *,float> distances;
        Cluster(BranchWithData *bwd)
        {
            branch = bwd;
            size = 1;
        }
        Cluster(Cluster *_U, Cluster *_V)
        {
            U = _U;
            V = _V;
            size = U->size + V->size;
        }
        void to_branch_data(std::vector<Branch *> &branches)
        {
            if (branch)
                branches.push_back(branch->b);
            else if (U && V)
            {
                U->to_branch_data(branches);
                V->to_branch_data(branches);
            }
        }
        float ward_dist(Cluster *B)
        {
            auto it = distances.find(B);
            if (it!=distances.end())
                return it->second;
            it = B->distances.find(this);
            if (it != B->distances.end())
                return it->second;
            if (branch)
            {//dist between single branch and cluster. Reduce second cluster
                if (B->branch)
                {//dist between single branches.
                    float distance = Clusterizer::dist(*branch, *(B->branch)).from;
                    distances.emplace(B, distance);
                    return distance;
                } 
                else
                {
                    float d1 = B->U->ward_dist(this);
                    float d2 = B->V->ward_dist(this);
                    float d3 = B->U->ward_dist(B->V);
                    float a = (float)(size + B->U->size)/(size + B->size);
                    float b = (float)(size + B->V->size)/(size + B->size);
                    float c = -(float)size/(size + B->size);
                    float distance = a*d1 + b*d2 + c*d3;
                    distances.emplace(B, distance);
                    return distance;
                }
            }
            else
            {//dist between two clusters. Reduce first cluster
                float d1 = U->ward_dist(B);
                float d2 = V->ward_dist(B);
                float d3 = U->ward_dist(V);
                float a = (float)(B->size + U->size)/(size + B->size);
                float b = (float)(B->size + V->size)/(size + B->size);
                float c = -(float)(B->size)/(size + B->size);
                float distance = a*d1 + b*d2 + c*d3;
                distances.emplace(B, distance);
                //fprintf(stderr,"dist calc %f*%f %f*%f %f*%f \n",d1,a,d2,b,d3,c);
                return distance;
            }
        }
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
    static bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, float min, float max);
    static bool match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc, float min, float max);
    static Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    static Answer dist_simple(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    static Answer dist_slow(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    static Answer dist_Nsection(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    Answer cluster_dist_min(Cluster &c1, Cluster &c2, float min = 1.0, float max = 0.0);
    Clusterizer()
    {
    }
    static float delta;
    BranchHeap branchHeap;
    std::vector<BranchWithData> branches;
    DistDataTable ddt;
    static std::vector<float> weights;
}; 