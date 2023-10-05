#pragma once
#include "clustering.h"
#include "hierarcial_clustering.h"

class DummyClusteringHelper : public ClusteringHelper
{
public:
    virtual BranchClusteringData *convert_branch(Block &settings, Branch *base, ClusteringContext *ctx, BaseBranchClusteringData &data) override
    {
        return new BranchClusteringData();
    }
    virtual void clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx) override
    {
        delete base;
    }
    virtual IntermediateClusteringData *prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                  ClusteringContext *ctx) override
    {
        IntermediateClusteringDataDDT *data = new IntermediateClusteringDataDDT();
        data->ddt.create(branches.size());
        for (int i=0;i<branches.size();i++)
        {
            for (int j=0;j<branches.size();j++)
            {
                float d = urand(0,3);
                if (d > 1.5)
                    d = 1e9;
                Answer a = Answer(true,d,d);
                DistData dd = DistData(d,0);
                data->ddt.set(i,j,a,dd);
            }
        }
        return data;
    }
};

class DummyClusteringBase : public ClusteringBase
{
public:
    virtual bool clusterize(Block &settings, IntermediateClusteringData *data, std::vector<ClusterStruct> &result) override
    {
        for (int i=0;i<data->elements_count;i+=2)
        {
            result.emplace_back();
            result.back().center = i;
            result.back().members = {std::pair<int, Transform>(i,Transform())};
            if (i+1<data->elements_count)
            {
                result.back().members.push_back(std::pair<int, Transform>(i+1,Transform()));
            }
        }
        return true;
    }
};