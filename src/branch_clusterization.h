#pragma once
#include "tree.h"
#include <vector>
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
        Answer *data;
        int n;
        void create(int _n)
        {
            n = _n;
            data = new Answer[n*n];
        }
        void clear()
        {
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
        std::vector<BranchWithData> &branches;
    };
    struct 
    bool set_branches(Tree &t, int layer);
    bool match_joints(Branch *b1, Branch *b2, std::vector<float> &matches, std::vector<int> &jc, float min, float max);
    bool match_child_branches(Joint *j1, Joint *j2, std::vector<float> &matches, std::vector<int> &jc, float min, float max);
    Answer dist(BranchWithData &bwd1, BranchWithData &bwd2, float min = 1.0, float max = 0.0);
    Answer cluster_dist_min(Cluster &c1, Cluster &c2, float min = 1.0, float max = 0.0);
    void 
    float delta = 0.5;
    BranchHeap branchHeap;
    std::vector<BranchWithData> branches;
    DistDataTable ddt;
    std::vector<float> weights{0.6, 0.25, 0.1, 0.035, 0.015};
}; 