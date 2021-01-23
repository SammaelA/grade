#pragma once
#include <list>
#include <vector>
#include "tinyEngine/utility/model.h"
#include "texture_atlas.h"
#include <glm/glm.hpp>
#include "tinyEngine/bbox.h"
    struct Billboard
    {
        int id = -1;
        int branch_id = -1;
        std::vector<glm::vec3> positions;
        glm::vec4 planeCoef; //billboard is always a plane ax+by+cz+d = 0 len(a,b,c) = 1
        bool instancing;
        void to_model(Model *m, TextureAtlas &atlas);
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
struct BillboardCloudData
{
    struct BillboardData
    {
        Billboard billboard;
        std::vector<glm::mat4> transforms; 
    };
    bool valid = false;
    std::vector<BillboardData> billboards;
    TextureAtlas atlas;
    BillboardCloudData() {};
};