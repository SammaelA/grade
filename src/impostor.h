#pragma once
#include "billboard_cloud.h"
#include "tree.h"


class ImpostorBaker : public BillboardCloudRaw
{
public:
    struct ImpostorGenerationParams
    {
        Quality quality;
        int slices_n = 8;
        bool need_top_view = true;
        bool fixed_colors = false;
        float add_rotation_y = 0;
        float leaf_size_mult = 1.15;
        int level_from = 0;
        int level_to = 1000;
    }; 
    ImpostorBaker(Quality _quality, int branch_level, std::vector<ClusterData> &clusters,
                                     std::vector<TreeTypeData> &_ttd, ImpostorsData *data)
    {
        ImpostorGenerationParams params;
        params.quality = _quality;
        quality = _quality;
        ttd = _ttd;
        prepare(branch_level,clusters,params,data);
    }
    ImpostorBaker(int tex_w, int tex_h, std::vector<TreeTypeData> &ttd) : BillboardCloudRaw(tex_w,tex_h,ttd) {};
    ImpostorBaker() {};
    void prepare(int branch_level, std::vector<ClusterData> &clusters,
                 ImpostorGenerationParams &params, ImpostorsData *data = nullptr);
    void prepare(ImpostorGenerationParams params, int branch_level, ClusterData &cluster, std::vector<TreeTypeData> &_ttd,
                 ImpostorsData *data, std::list<Impostor>::iterator &impostor);
    void prepare(Quality quality, int branch_level, ClusterData &cluster, std::vector<TreeTypeData> &_ttd,
                 ImpostorsData *data, std::list<Impostor>::iterator &impostor);
    void prepare_all_grove(GroveGenerationData &ggd, int branch_level, std::vector<ClusterData> &clusters,
                           ImpostorsData *data = nullptr);
    void make_impostor(Branch &b, Impostor &imp, ImpostorGenerationParams &params);
private:

};
class ImpostorRenderer
{
    public:
    friend class GroveRenderer;
    ImpostorRenderer(ImpostorsData *data = nullptr);
    ~ImpostorRenderer();
    void render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light, 
                glm::mat4 &shadow_tr, GLuint shadow_tex,
                glm::vec3 camera_pos = glm::vec3(0,0,0),
                glm::vec4 screen_size = glm::vec4(800,600,1/800,1/600), 
                bool to_shadow = false,
                GroveRendererDebugParams dbgpar = GroveRendererDebugParams());
    private:
    GLuint slicesBuffer = 0;
    GLuint impostorsDataBuffer = 0;
    ImpostorsData *data= nullptr;

    Shader impostorRenderer;
    Shader impostorRendererInstancing;
    
    std::vector<Model *> models;
    std::vector<uint> offsets;

    float hth = 0,vth = 0;
};