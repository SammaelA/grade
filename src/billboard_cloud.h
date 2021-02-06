#pragma once
#include "texture_atlas.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/instance.h"
#include "tree.h"
#include "branch_clusterization.h"
#include "billboard_cloud_data.h"
#include <list>
#include <vector>
class Visualizer;
class BillboardCloudRaw
{

public:
    enum RenderMode
    {
        NOTHING,
        ONLY_SINGLE,
        ONLY_INSTANCES,
        BOTH
    };
    BillboardCloudRaw(int tex_w, int tex_h, std::vector<TreeTypeData> &ttd);
    ~BillboardCloudRaw();
    void setup_preparation();
    void prepare(Tree &t, int branch_level, std::vector<Branch> &branches);
    void prepare(Tree &t, int branch_level, int layer);
    void prepare(Tree &t, int branch_level, std::vector<Clusterizer::Cluster> &clusters, std::list<int> &numbers, BillboardCloudData *data = nullptr);
    void render(glm::mat4 &projectionCamera);
    void set_textures(Texture _wood);
    static BBox get_bbox(Branch *branch, glm::vec3 a, glm::vec3 b, glm::vec3 c);
    void set_render_mode(RenderMode m)
    {
        renderMode = m;
    }
private:
    struct BranchProjectionData
    {
        float projection_err = 0.0;
        Branch *br;
        int parent_n;
        BranchProjectionData(float err, Branch *b, int p_n)
        {
            projection_err = err;
            br = b;
            parent_n = p_n;
        }
    };
    struct BillboardBox
    {
        Branch *b;
        BBox min_bbox;
        glm::vec3 base_joint;
        int parent;
        BillboardBox(Branch *_b, BBox &_bbox, glm::vec3 _base_joint, int par = -1) : min_bbox(_bbox)
        {
            b = _b;
            base_joint = _base_joint;
            parent = par;
        }
    };
    void prepare_branch(Tree &t, Branch *b, BBox &min_box, Visualizer &tg, int billboards_count);
    void create_billboard(TreeTypeData &ttd, Branch *b, BBox &min_box, Visualizer &tg, int id, Billboard &bill);
    BBox get_minimal_bbox(Branch *b);
    static void update_bbox(Branch *branch, glm::mat4 &rot, glm::vec4 &mn, glm::vec4 &mx);
    glm::mat4 get_viewproj(BBox &b);
    static bool BPD_comp(BranchProjectionData &a, BranchProjectionData &b);
    float projection_error_rec(Branch *b, glm::vec3 &n, float d);
    int billboard_count = 256;
    bool ready = false;
    TextureAtlas atlas;
    std::vector<TreeTypeData> ttd;
    Shader rendererToTexture;
    Shader billboardRenderer;
    Shader billboardRendererInstancing;
    Model *cloud;
    std::vector<Instance *> instances;
    Texture wood;
    std::vector<Billboard> billboards;
    RenderMode renderMode = ONLY_SINGLE;
};
struct BillboardCloudRenderer
{   
    enum RenderMode
    {
        NOTHING,
        ONLY_SINGLE,
        ONLY_INSTANCES,
        BOTH
    };
    BillboardCloudRenderer(BillboardCloudData *data = nullptr);
    ~BillboardCloudRenderer();
    void render(glm::mat4 &projectionCamera, glm::vec3 camera_pos = glm::vec3(0,0,0), glm::vec2 LOD_min_max = glm::vec2(0,0),
                glm::vec4 screen_size = glm::vec4(800,600,1/800,1/600));
    void set_render_mode(RenderMode m)
    {
        renderMode = m;
    }
    private:
    BillboardCloudData *data= nullptr;
    Shader rendererToTexture;
    Shader billboardRenderer;
    Shader billboardRendererInstancing;
    Model *cloud = nullptr;
    std::vector<Instance *> instances;
    RenderMode renderMode = ONLY_SINGLE;

};