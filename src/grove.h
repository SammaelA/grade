#pragma once
#include "tinyEngine/utility.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/instance.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility.h"
#include "billboard_cloud_data.h"
#include "timestamp.h"
#include "rendering_SSBO_structs.h"
#include "tinyEngine/camera.h"
#include "limits.h"
#include "tinyEngine/utility/bit_vector.h"
class ImpostorRenderer;
class BillboardCloudRenderer;
class GroveGenerationData;
struct BBox;
struct PackedLeaf
{
    std::vector<glm::vec3> edges;
};
struct PackedJoint
{
    glm::vec3 pos;
    float r;
    PackedJoint() { r = 0; }
    PackedJoint(glm::vec3 &_pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
};
struct PackedBranch
{
    std::vector<PackedJoint> joints;
    std::vector<PackedLeaf> leaves;
    std::vector<std::vector<float>> r_mults;
    int level;
    uint type_id;
    glm::vec4 plane_coef;
};
struct BranchCatalogue
{
    static const unsigned LEVEL_BITS = 4;
    std::vector<std::vector<PackedBranch>> branches;
    BranchCatalogue(int levels)
    {
        if (levels >= (1 << LEVEL_BITS))
        {
            logerr("Branch catalogue created with too many branch levels: %d. Max value is %d", levels, (1 << LEVEL_BITS) - 1);
            levels = (1 << LEVEL_BITS) - 1;
        }
        for (int i = 0; i < levels; i++)
        {
            branches.push_back(std::vector<PackedBranch>());
            occupied.push_back(BitVector());
            first_free.push_back(0);
        }
    }
    int levels() const
    {
        return branches.size();
    }
    PackedBranch &get(unsigned pos)
    {
        if (occupied[pos & ((1 << LEVEL_BITS) - 1)].get(pos >> LEVEL_BITS))
            return branches[pos & ((1 << LEVEL_BITS) - 1)][pos >> LEVEL_BITS];
        else
            return emptyBranch;
    }
    int add(PackedBranch &b, int level)
    {
        if (level < 0 || level >= branches.size())
            return -1;
        int pos = first_free[level];
        if (pos == branches[level].size())
        {
            //increase vector size
            first_free[level]++;
            occupied[level].push_back(true);
            branches[level].push_back(b);
        }
        else
        {
            branches[level][pos] = b;
            occupied[level].set_unsafe(pos,true);
            first_free[level] = branches[level].size(); 
            for (int i=pos;i<branches[level].size();i++)
            {
                if (occupied[level].get_unsafe(i))
                {
                    first_free[level] = i;
                    break;
                }
            }
        }
        return ((pos) << LEVEL_BITS) + level;
    }
    int get_level_size(int level)
    {
        if (level < 0)
            level = 0;
        if (level >= branches.size())
            level = branches.size() - 1;
        int sz = 0;
        for (int i=0;i<occupied[level].size();i++)
        {
            sz += occupied[level].get_unsafe(i);
        }
        return sz;
    }
    void remove(unsigned id)
    {
        int level = id & ((1 << LEVEL_BITS) - 1);
        int pos = id >> LEVEL_BITS;
        occupied[level].set(pos, false);
        first_free[level] = MIN(first_free[level],pos);
    }
    private:
    PackedBranch emptyBranch;
    std::vector<BitVector> occupied;
    std::vector<int> first_free;

};
struct BranchStructure
{
    unsigned pos;
    std::vector<BranchStructure> childBranches;
    std::vector<std::pair<glm::mat4, BranchStructure>> childBranchesInstanced;
    explicit BranchStructure(unsigned _pos = 0) {pos = _pos;}
};
struct InstancedBranch
{
    std::vector<unsigned> branches;
    InstanceDataArrays IDA;
    BBox bbox;
};
struct GrovePacked
{
    glm::vec3 center;
    std::string ggd_name;
    BranchCatalogue instancedCatalogue;

    std::list<InstancedBranch> instancedBranches;
    std::vector<BillboardCloudData> clouds;
    std::vector<ImpostorsData> impostors;
    GrovePacked() : instancedCatalogue(MAX_BRANCH_LEVELS){};
};
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
    struct Instance2 : InstancingReadyModel
    {
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
    GroveRenderer(GrovePacked *_source, GroveGenerationData *_ggd, int LODs_count, std::vector<float> &max_distances,
                  bool print_perf, Precision precision);
    GroveRenderer();
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
    GroveGenerationData *ggd;
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