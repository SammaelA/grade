#pragma once
#include "common_utils/utility.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
#include "common_utils/utility.h"
#include "core/billboard_cloud_data.h"
#include "timestamp.h"
#include "render/rendering_SSBO_structs.h"
#include "tinyEngine/camera.h"
#include "limits.h"
#include "common_utils/bit_vector.h"
#include "core/limits.h"
#include "core/grove.h"

class ImpostorRenderer;
struct BillboardCloudRenderer;

struct GroveRendererDebugParams
{
    bool need_focus_model = false;
    int model_focused = 0;
};
class GroveRenderer
{
public:
    static const int base_level = 1;
    enum Precision
    {
        LOW = 0,
        MEDIUM = 1,
        DEBUG = 2
    };
    struct Instance2
    {
        int id = -1;
        Model *m;
        DrawElementsIndirectCommand cmd;
        uint type;
        InstanceDataArrays &ida;
        BBox bbox;
        Instance2(Model *m, uint type, InstanceDataArrays &ida);
        ~Instance2();
    };
    struct LOD
    {
        std::vector<std::pair<uint,Model *>> models;
        BillboardCloudRenderer *cloud = nullptr;
        ImpostorRenderer *imp_rend = nullptr;
        std::vector<Instance2> instances;
        std::vector<Instance2> leaves_instances;
        float max_dist;
    };
    int get_max_LOD() {return MIN(4,LODs.size() - 1);}
    void render(int lod, glm::mat4 &projection, glm::mat4 &view, Camera &camera, glm::vec2 screen_size, DirectedLight &light, 
                GroveRendererDebugParams dbgpar, glm::mat4 &shadow_tr, GLuint shadow_tex, bool to_shadow = false);
    GroveRenderer(const GrovePacked *_source, AABB2D scene_bbox, const std::vector<TreeTypeData> &types, 
                  int LODs_count, std::vector<float> &max_distances, bool print_perf, Precision precision);
    ~GroveRenderer();
private:
    struct TypeDescriptionForRender
    {
        ImpostorRenderer *imp = nullptr;
        BillboardCloudRenderer *bill = nullptr;
        MultiDrawRendDesc rendDesc;
        
    };
    struct CellsParams
    {
        int x_cells;
        int y_cells;
        float x_size;
        float y_size;
        glm::vec3 start_pos;
    };
    void add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch, int up_to_level, bool need_leaves = false);
    void IDA_to_bufer(InstanceDataArrays &ida, std::vector<LodData> &lods, std::vector<InstanceData> &instances,
                      std::vector<ModelData> &models, std::vector<TypeData> &types, bool is_leaf = false);
    DrawElementsIndirectCommand model_to_base(Model *m, BBox &bb);
    void pack_bb_to_model_data(ModelData &md, BBox &bb);
    void prepare_wood_types_atlas();
    std::vector<LOD> LODs;
    std::vector<TypeDescriptionForRender> types_descs;
    Shader renderer;
    Shader rendererInstancing;
    Shader lodCompute;
    Shader cellsCompute;
    Shader shadowRendererInstancing;
    GrovePacked *source;
    AABB2D scene_bbox;
    std::vector<TreeTypeData> types;
    GLuint lods_buf, instances_buf, models_buf, types_buf;
    GLuint cur_insts_buf, cur_models_buf, cur_types_buf;
    GLuint draw_indirect_buffer, cells_buf;
    int total_models_count = 0;
    Timestamp ts;
    Model *base_container;
    TextureAtlas *atlas = nullptr;
    Camera prev_camera;
    uint64_t frames = 0;
    CellsParams cellsInfo;

};