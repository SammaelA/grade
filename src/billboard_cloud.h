#pragma once
#include "texture_atlas.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility/model.h"
#include "tree.h"
#include <list>
#include <vector>
class TreeGenerator;
class Visualizer;
class BillboardCloud
{

    public:
    struct BBox
    {
        glm::vec3 position;
        glm::vec3 sizes;
        glm::vec3 a;
        glm::vec3 b;
        glm::vec3 c;
        float V(){return sizes.x*sizes.y*sizes.z;}
        bool special = false;
    };
    BillboardCloud(int tex_w, int tex_h);
    ~BillboardCloud();
    void setup_preparation();
    void prepare(Tree &t, int layer);
    void render(glm::mat4 &projectionCamera);
    void set_textures(Texture *wood);
    static BBox get_bbox(Branch *branch, glm::vec3 a, glm::vec3 b, glm::vec3 c);

    private:
    struct Billboard
    {
        int id = - 1;
        std::vector<glm::vec3> positions;
        glm::vec4 planeCoef; //billboard is always a plane ax+by+cz+d = 0 len(a,b,c) = 1
        void to_model(Model *m, TextureAtlas &atlas);
        Billboard() {};
        Billboard(const Billboard &b)
        {
            this->id = b.id;
            this->positions = b.positions;
        }
        Billboard(const BBox &box, int id, int type, glm::vec3 base_joint);
    };
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
        BillboardBox(Branch *_b, BBox &_bbox, glm::vec3 _base_joint, int par = -1):
        min_bbox(_bbox)
        {
            b = _b;
            base_joint = _base_joint;
            parent = par;
        }
    };
    void prepare_branch(Tree &t, Branch *b, BBox &min_box, Visualizer &tg, int billboards_count);
    void create_billboard(Tree &t, Branch *b, BBox &min_box, Visualizer &tg, int id, Billboard &bill);
    BBox get_minimal_bbox(Branch *b);
    static void update_bbox(Branch *branch, glm::mat4 &rot, glm::vec4 &mn, glm::vec4 &mx);
    glm::mat4 get_viewproj(BBox &b);
    static bool BPD_comp(BranchProjectionData &a, BranchProjectionData &b);
    float projection_error_rec(Branch *b, glm::vec3 &n, float d);

    int billboard_count = 256;
    bool ready = false;
    TextureAtlas atlas;
    Shader rendererToTexture;
    Shader billboardRenderer;
    Model *cloud;
    Texture *pwood = nullptr;
    std::vector<Billboard> billboards;
};