#pragma once
#include <list>
#include <vector>
#include "tinyEngine/model.h"
#include "graphics_utils/texture_atlas.h"
#include <glm/glm.hpp>
#include "common_utils/bbox.h"
#include "save_utils/serialization.h"

struct InstanceDataArrays
{
    friend class boost::serialization::access;

    std::vector<glm::mat4> transforms;
    std::vector<glm::vec3> centers_par;
    std::vector<glm::vec3> centers_self;
    std::vector<int> type_ids;
    std::vector<int> tree_ids;

    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & transforms;
      ar & centers_par;
      ar & centers_self;
      ar & type_ids;
      ar & tree_ids;
    }
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
        friend class boost::serialization::access;

        int id = -1;
        int branch_id = -1;
        std::vector<glm::vec3> positions;
        glm::vec4 planeCoef; //billboard is always a plane ax+by+cz+d = 0 len(a,b,c) = 1
        bool instancing;
        void to_model(Model *m, const TextureAtlas &atlas);
        std::vector<glm::vec3> get_tc(const TextureAtlas &atlas);
        Billboard(){};
        Billboard(const Billboard &b)
        {
            this->id = b.id;
            this->branch_id = b.branch_id;
            this->positions = b.positions;
            this->instancing = b.instancing;
        }
        Billboard(const BBox &box, int id, int branch_id, int type, glm::vec3 base_joint, bool _instancing = false);
      
    private:
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar & id;
        ar & branch_id;
        ar & positions;
        ar & planeCoef;
        ar & instancing;
      }
    };
struct BillboardData
{
    glm::vec3 base_position;
    std::vector<Billboard> billboards;
    InstanceDataArrays IDA;
    int id = -1;
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
    friend class boost::serialization::access;

    glm::vec3 center;
    float r;
    float h_2;//distance from center to base, h/2;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & center;
      ar & r;
      ar & h_2;
    }
};
struct Impostor
{
    friend class boost::serialization::access;

    int id = -1;
    std::vector<Billboard> slices;
    Billboard top_slice;
    BCyl bcyl;
    InstanceDataArrays IDA;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & id;
      ar & slices;
      ar & top_slice;
      ar & bcyl;
      ar & IDA;
    }
};
struct ImpostorsData
{ 
    friend class boost::serialization::access;

    bool valid = false;
    int level = 0;
    std::list<Impostor> impostors;
    TextureAtlas atlas;
    ImpostorsData() {};
    ~ImpostorsData() {atlas.destroy();}
  
  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & valid;
      ar & level;
      ar & impostors;
      ar & atlas;
    }
};