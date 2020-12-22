#pragma once
class Mark
{
    public:
    #define UNKNOWN 0
    int type() {return UNKNOWN;}
    virtual ~Mark() = 0;
};
class ClusteringJointMark : Mark
{
    public:
    #define CLUSTERING 1
    int type() {return CLUSTERING;}
    int joints_count;
    virtual ~ClusteringJointMark() {};
};