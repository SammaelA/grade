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
    struct Cluster
    {
        std::vector<Branch *> branches;
    };
    bool set_branches(Tree &t, int layer);
    Answer dist(Branch *b1, Branch *b2, float min = 1.0, float max = 0.0);
    float delta = 0.5;
    BranchHeap branchHeap;
}; 