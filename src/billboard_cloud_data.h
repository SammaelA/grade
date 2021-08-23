#pragma once
#include <list>
#include <vector>
#include "tinyEngine/utility/model.h"
#include "texture_atlas.h"
#include <glm/glm.hpp>
#include "tinyEngine/bbox.h"
struct InstanceDataArrays
{
    std::vector<glm::mat4> transforms;
    std::vector<glm::vec3> centers_par;
    std::vector<glm::vec3> centers_self;
    std::vector<int> type_ids;
    std::vector<int> tree_ids;
};
struct InstancingReadyModel
{
    int id = -1;
};
    struct MultiDrawRendDesc
    {
        int type_id;
        int cmd_buffer_offset;
        int current_types_offset;
        int max_models;
        int cmd_size;
        int base_vertex_id;
    };
    struct Billboard
    {
        int id = -1;
        int branch_id = -1;
        std::vector<glm::vec3> positions;
        glm::vec4 planeCoef; //billboard is always a plane ax+by+cz+d = 0 len(a,b,c) = 1
        bool instancing;
        void to_model(Model *m, TextureAtlas &atlas);
        std::vector<glm::vec3> get_tc(TextureAtlas &atlas);
        Billboard(){};
        Billboard(const Billboard &b)
        {
            this->id = b.id;
            this->branch_id = b.branch_id;
            this->positions = b.positions;
            this->instancing = b.instancing;
        }
        Billboard(const BBox &box, int id, int branch_id, int type, glm::vec3 base_joint, bool _instancing = false);
    };
struct BillboardData : InstancingReadyModel
{
    glm::vec3 base_position;
    std::vector<Billboard> billboards;
    InstanceDataArrays IDA;
};
struct BillboardCloudData
{
    bool valid = false;
    int level = 0;
    std::list<BillboardData> billboards;
    TextureAtlas atlas;
    BillboardCloudData() {};
    ~BillboardCloudData() {atlas.destroy();}
};
struct BCyl
{
    glm::vec3 center;
    float r;
    float h_2;//distance from center to base, h/2;
};
struct Impostor : InstancingReadyModel
{
    std::vector<Billboard> slices;
    Billboard top_slice;
    BCyl bcyl;
    InstanceDataArrays IDA;
};
struct ImpostorsData
{
    bool valid = false;
    int level = 0;
    std::vector<Impostor> impostors;
    TextureAtlas atlas;
    ImpostorsData() {};
    ~ImpostorsData() {atlas.destroy();}
};