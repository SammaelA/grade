#pragma once
#include "graphics_utils/billboard_cloud.h"
#include "core/tree.h"

class ImpostorBaker : public BillboardCloudRaw
{
public:
    struct ImpostorGenerationParams : public BillboardCloudRaw::BillboardGenerationParams
    {
        int slices_n = 8;
        bool need_top_view = true;
        float add_rotation_y = 0;
    }; 
    ImpostorBaker() {};
    void prepare(int branch_level, std::vector<ClusterData> &clusters,
                 ImpostorGenerationParams &params, ImpostorsData *data = nullptr);
    void prepare(ImpostorGenerationParams params, int branch_level, ClusterData &cluster, const std::vector<TreeTypeData> &_ttd,
                 ImpostorsData *data, std::list<Impostor>::iterator &impostor, int clusters_expected = 1);
    void prepare(Quality quality, int branch_level, ClusterData &cluster, const std::vector<TreeTypeData> &_ttd,
                 ImpostorsData *data, std::list<Impostor>::iterator &impostor);
    void make_impostor(Branch &b, const TreeTypeData &tree_type, Impostor &imp, ImpostorGenerationParams &params, 
                       TextureAtlas &atlas, BBox &bbox);
private:

};